#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
artifact="${1:-$repo_root/vendor/rsutils/wheels-linux-$(uname -m).zip}"
if [ ! -f "$artifact" ]; then
    printf 'RS-utils package not found: %s\nUse a supported bundled wheel archive, or pass an instructor ZIP/wheel path.\n' "$artifact" >&2
    exit 1
fi
artifact="$(realpath -- "$artifact")"
wheel_bundle=false
case "$artifact" in
    *.zip)
        actual_sha="$(sha256sum "$artifact" | cut -d ' ' -f 1)"
        case "$actual_sha" in
            f071a936745d1dfc30e597ecf14a86aff2e3788db746e494fb5dc0dc02c75b80) ;;
            63f1c9b13462c43127d0c375ddda6c1f0749cbd2df96669b7260576250efc826|e96353595e8aa767b947e439d89a105351bb45703c0e3974009f64d385b3492b) wheel_bundle=true ;;
            *) printf 'This ZIP does not match an instructor release-prep source or wheel bundle. Obtain the matching package or update the pinned checksum.\n' >&2; exit 1 ;;
        esac
        ;;
    *.whl) ;;
    *) printf 'Pass the release-prep source ZIP, wheel-bundle ZIP, or a compatible Linux wheel.\n' >&2; exit 1 ;;
esac

python3 -c 'import sys; assert sys.version_info >= (3, 11), "RS-utils needs Python 3.11 or newer; rebuild the Ubuntu 24.04 container"'
venv="$repo_root/.venv-rsutils"
python3 -m venv --system-site-packages "$venv"
"$venv/bin/python" -m pip install --upgrade pip packaging
if [ "$wheel_bundle" = true ]; then
    wheel_dir="$(mktemp -d)"
    trap 'rm -rf "$wheel_dir"' EXIT
    artifact="$("$venv/bin/python" - "$artifact" "$wheel_dir" <<'PYWHEEL'
from pathlib import Path
from packaging.tags import sys_tags
from packaging.utils import parse_wheel_filename
import sys, zipfile
ranks = {tag: rank for rank, tag in enumerate(sys_tags())}
choices = []
with zipfile.ZipFile(sys.argv[1]) as bundle:
    for member in bundle.namelist():
        if not member.endswith('.whl'):
            continue
        name, version, _, tags = parse_wheel_filename(Path(member).name)
        compatible = tags & ranks.keys()
        if name == 'rslmto-utils' and str(version) == '0.11.0rc1' and compatible:
            choices.append((min(ranks[tag] for tag in compatible), member))
    if not choices:
        raise SystemExit('No wheel in this bundle matches this Python/CPU/platform. Use the matching architecture bundle.')
    member = min(choices)[1]
    target = Path(sys.argv[2]) / Path(member).name
    target.write_bytes(bundle.read(member))
print(target)
PYWHEEL
)"
    printf 'Selected %s\n' "$(basename "$artifact")"
fi
# Remove both distribution names so even an equal-version snapshot is replaced.
"$venv/bin/python" -m pip uninstall --yes rsutil rslmto-utils
"$venv/bin/python" -m pip install "$artifact"
"$venv/bin/python" - <<'PY'
from importlib.metadata import version
import rsutil
installed = version('rslmto-utils')
if installed != '0.11.0rc1':
    raise RuntimeError(f'Expected instructor release-prep 0.11.0rc1, found {installed}')
if not rsutil.native_library_available():
    raise RuntimeError('RS-utils installed, but its native Fortran library could not be loaded')
print(f'rslmto-utils {installed}: Python import and native Fortran library ready')
PY
"$venv/bin/rsutil" --help >/dev/null
printf '\nReady. Open a new terminal, or run:\n  source "%s/bin/activate"\n  rsutil --help\n' "$venv"
