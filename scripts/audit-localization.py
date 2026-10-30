#!/usr/bin/env python3
"""Audit UI literal translation call sites (not a completeness proof)."""
import pathlib
import re
import sys
root=pathlib.Path(__file__).resolve().parents[1]
pattern=re.compile(r'(?:\bui\.(?:label|button|menu_button|strong|small|heading)\(\s*|(?:RichText|Button|Label)::new\(\s*|\.(?:on_hover_text|hint_text)\(\s*)"([^"\n]*[A-Za-z][^"\n]*)"')
count=0
allowed={'AppimageGene','Language / 语言：','°C'}
unreviewed=0
for path in sorted((root/'crates/hx-app/src').rglob('*.rs')):
    if path.name=='i18n.rs':continue
    text=path.read_text()
    for match in pattern.finditer(text):
        count+=1
        known=match[1] in allowed
        if not known:unreviewed+=1
        print(f'{path.relative_to(root)}:{text[:match.start()].count(chr(10))+1}: {match[1]}'+(' [allowed: identity/language name/unit]' if known else ' [REVIEW]'))
print(f'{unreviewed} unreviewed direct literal sites; {count} total. Dynamic/indirect text needs separate review.')
if '--check' in sys.argv and unreviewed:raise SystemExit(1)
