#!/usr/bin/env python3
"""Derive a stable positive 31-bit seed from biological identifiers.

The scheme deliberately excludes scheduler task numbers, manifest row numbers,
batch sizes, and filesystem paths.  Its byte-level definition is frozen as
``stable-seed-v1`` so that independent implementations can reproduce it.
"""

from __future__ import annotations

import argparse
import hashlib
import json


SCHEME = "stable-seed-v1"
MAX_SEED = 2_147_483_646


def stable_seed(base: int, namespace: str, components: list[str]) -> int:
    payload = json.dumps(
        [SCHEME, base, namespace, *components],
        ensure_ascii=False,
        separators=(",", ":"),
    ).encode("utf-8")
    digest = hashlib.sha256(payload).digest()
    return 1 + (int.from_bytes(digest[:8], "big") % MAX_SEED)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=int, default=13)
    parser.add_argument("--namespace", required=True)
    parser.add_argument("components", nargs="+")
    args = parser.parse_args()
    print(stable_seed(args.base, args.namespace, args.components))


if __name__ == "__main__":
    main()
