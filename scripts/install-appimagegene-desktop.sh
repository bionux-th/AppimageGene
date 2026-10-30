#!/usr/bin/env bash
# Product desktop integration installer.
# Install current user's launcher/icon only. Does not claim MIME defaults.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
data=${XDG_DATA_HOME:-$HOME/.local/share}
mkdir -p "$data/applications" "$data/icons/hicolor/scalable/apps" "$data/icons/hicolor/256x256/apps"
install -m 644 "$root/assets/icons/appimagegene.svg" "$data/icons/hicolor/scalable/apps/appimagegene.svg"
install -m 644 "$root/assets/icons/appimagegene.png" "$data/icons/hicolor/256x256/apps/appimagegene.png"
python3 - "$root" "$data/applications/appimagegene.desktop" <<'PY'
import pathlib,sys
root=pathlib.Path(sys.argv[1]); destination=pathlib.Path(sys.argv[2])
binary=str(root/'target/release/appimagegene')
# Escape Desktop Entry quoted Exec arguments, including reserved characters.
escaped=binary.replace('\\','\\\\').replace('"','\\"').replace('`','\\`').replace('$','\\$').replace('%','%%')
text=(root/'assets/appimagegene.desktop').read_text()
text=text.replace('Exec=appimagegene %F',f'Exec="{escaped}" %F')
destination.write_text(text)
print('Installed launcher:',destination)
PY
if command -v update-desktop-database >/dev/null; then update-desktop-database "$data/applications"; fi
if command -v gtk-update-icon-cache >/dev/null && [ -f "$data/icons/hicolor/index.theme" ]; then gtk-update-icon-cache -f -t "$data/icons/hicolor"; fi
if command -v kbuildsycoca6 >/dev/null; then kbuildsycoca6 --noincremental; fi
