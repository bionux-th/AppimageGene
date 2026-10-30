#!/usr/bin/env python3
"""Package only a built executable and its actual MinGW DLL imports; no downloads."""
import argparse
import hashlib
import pathlib
import shutil
import subprocess
import zipfile

parser=argparse.ArgumentParser()
parser.add_argument('--target',default='x86_64-pc-windows-msvc')
args=parser.parse_args()
root=pathlib.Path(__file__).resolve().parents[1]
binary=root/'target'/args.target/'release/appimagegene.exe'
if not binary.is_file():raise SystemExit('Build appimagegene.exe first: '+str(binary))
out=root/'dist'/('AppimageGene-windows-x64-'+args.target.rsplit('-',1)[-1])
out.mkdir(parents=True,exist_ok=True)
shutil.copy2(binary,out/binary.name)
shutil.copytree(root/'assets/resources',out/'assets/resources',dirs_exist_ok=True)
if args.target.endswith('-gnu'):
    pending=[out/binary.name];seen=set()
    while pending:
        image=pending.pop()
        text=subprocess.check_output(['x86_64-w64-mingw32-objdump','-p',str(image)],text=True)
        for line in text.splitlines():
            if 'DLL Name:' not in line:continue
            name=line.split('DLL Name:',1)[1].strip()
            if name.lower()in seen:continue
            seen.add(name.lower())
            location=subprocess.check_output(['x86_64-w64-mingw32-gcc','-print-file-name='+name],text=True).strip()
            source=pathlib.Path(location)
            if not source.is_file():source=pathlib.Path('/usr/x86_64-w64-mingw32/bin')/name
            if source.is_file():shutil.copy2(source,out/name);pending.append(out/name)
shutil.copy2(root/'README.md',out/'README.md')
licenses=out/'licenses';licenses.mkdir(exist_ok=True)
shutil.copy2(root/'vendor/viennarna/COPYING',licenses/'ViennaRNA-COPYING')
shutil.copy2(root/'vendor/viennarna/AUTHORS',licenses/'ViennaRNA-AUTHORS')
if args.target.endswith('-gnu'):
    for path in pathlib.Path('/usr/share/licenses/mingw-w64-crt').glob('COPYING*'):
        shutil.copy2(path,licenses/path.name)
for name in ['LICENSE','LICENSE-MIT','LICENSE-APACHE']:
    if(root/name).is_file():shutil.copy2(root/name,out/name)
archive=out.with_suffix('.zip')
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED)as z:
    for path in sorted(out.rglob('*')):
        if path.is_file():z.write(path,path.relative_to(out.parent))
print(archive)
print('SHA256:',hashlib.sha256(archive.read_bytes()).hexdigest())
