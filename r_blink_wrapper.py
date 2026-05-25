"""Python wrapper to prepare inputs and call `run_blink.R`.
Converts PLINK binary to a PLINK --recode A raw table if needed, then calls Rscript.
"""
import os
import sys
import subprocess
import tempfile
from pathlib import Path

def find_plink():
    from shutil import which
    plink = which('plink1') or which('p-link') or which('plink2') or which('plink') or which('plink_mac')
    if plink:
        return plink
    cwd = Path.cwd()
    candidates = [
        cwd / 'plink',
        cwd / 'plink1',
        cwd / 'p-link',
        cwd / 'plink2',
        cwd / 'plink_mac',
        cwd / 'plink.exe',
    ]
    for cand in candidates:
        if cand.exists() and cand.is_file() and os.access(cand, os.X_OK):
            # Check if it's executable for the current OS/architecture
            try:
                # Run a quick plink --version or similar to ensure it's not a Mac binary in Linux
                subprocess.run([str(cand), '--noweb'], capture_output=True, check=False)
                return str(cand)
            except OSError:
                pass
    return None

def find_rscript(r_executable='Rscript'):
    from shutil import which
    return which(str(r_executable))

def plink_to_raw(plink_prefix, out_prefix=None):
    plink = find_plink()
    if plink is None:
        raise FileNotFoundError('PLINK executable not found in PATH or workspace')
    if out_prefix is None:
        out_prefix = plink_prefix + '_to_raw'
    raw_path = out_prefix + '.raw'
    # Reuse existing conversion when available to avoid repeating expensive PLINK work.
    if os.path.exists(raw_path):
        return raw_path
    cmd = [plink, '--noweb', '--bfile', plink_prefix, '--recodeA', '--out', out_prefix]
    subprocess.run(cmd, check=True)
    if not os.path.exists(raw_path):
        raise FileNotFoundError(f'PLINK conversion failed, expected {raw_path}')
    return raw_path

def call_r_blink(geno_csv, pheno_csv, covar_csv, out_prefix, trait=None, chunk_size=2000, method='blink', r_executable='Rscript', strict_native=True):
    rscript_path = find_rscript(r_executable)
    if not rscript_path:
        raise FileNotFoundError(f'{r_executable} not found in PATH; please install R and ensure it is available')
    script_path = Path(__file__).parent / 'run_blink.R'
    if not script_path.exists():
        # fallback to repository root location
        script_path = Path('run_blink.R')
    cmd = [
        rscript_path, str(script_path),
        '--geno_csv', str(geno_csv),
        '--pheno', str(pheno_csv),
        '--out', str(out_prefix),
        '--chunk_size', str(int(chunk_size)),
        '--method', str(method).lower(),
    ]
    if strict_native:
        cmd.extend(['--strict_native', '1'])
    if covar_csv:
        cmd.extend(['--covar', str(covar_csv)])
    if trait:
        cmd.extend(['--trait', str(trait)])
    print('Calling Rscript:', ' '.join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    print(proc.stdout)
    if proc.returncode != 0:
        print(proc.stderr, file=sys.stderr)
        raise RuntimeError(f'Rscript run failed with code {proc.returncode}')
    return proc

def run(geno_path, geno_format, pheno_csv, covar_csv, out_prefix, trait=None, chunk_size=2000, method='blink', r_executable='Rscript', strict_native=True):
    # geno_format: 'plink' or 'csv'
    geno_csv = geno_path
    rscript_path = find_rscript(r_executable)
    if not rscript_path:
        raise FileNotFoundError(f'{r_executable} not found in PATH; please install R and ensure it is available')
    if geno_format.lower() == 'plink':
        # Convert once and reuse the cached RAW across models for speed.
        cache_key = str(geno_path).replace(os.sep, '_').replace(':', '_')
        cache_prefix = str(Path(out_prefix).parent / f"_plink_raw_cache_{cache_key}")
        print('Converting PLINK to raw format (cached, may take time on first run)...')
        geno_csv = plink_to_raw(geno_path, cache_prefix)
    elif geno_format.lower() == 'csv':
        geno_csv = geno_path
    else:
        raise ValueError('Unsupported GENO format: ' + str(geno_format))
    return call_r_blink(
        geno_csv,
        pheno_csv,
        covar_csv,
        out_prefix,
        trait=trait,
        chunk_size=chunk_size,
        method=method,
        r_executable=r_executable,
        strict_native=strict_native,
    )

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--geno', required=True)
    p.add_argument('--format', default='plink')
    p.add_argument('--pheno', required=True)
    p.add_argument('--covar', default='')
    p.add_argument('--out', required=True)
    p.add_argument('--trait', default='')
    p.add_argument('--chunk', type=int, default=2000)
    p.add_argument('--method', default='blink')
    p.add_argument('--rscript', default='Rscript')
    p.add_argument('--strict-native', action='store_true', help='Require native package method; disallow linear fallback')
    p.add_argument('--allow-fallback', action='store_true', help='Allow linear fallback if native package fails/unavailable')
    args = p.parse_args()
    run(
        args.geno,
        args.format,
        args.pheno,
        args.covar or None,
        args.out,
        trait=args.trait or None,
        chunk_size=args.chunk,
        method=args.method,
        r_executable=args.rscript,
        strict_native=(args.strict_native or (not args.allow_fallback)),
    )
