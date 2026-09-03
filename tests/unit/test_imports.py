"""Every module imports, including the ones only a real run reaches."""

from __future__ import annotations

import importlib
import pkgutil

import pytest

import pymopsmap

MODULES = sorted(
    name
    for _, name, _ in pkgutil.walk_packages(
        pymopsmap.__path__, prefix="pymopsmap."
    )
)


@pytest.mark.parametrize("name", MODULES)
def test_the_module_imports(name):
    importlib.import_module(name)


def test_deferred_imports_resolve():
    """
    Functions import their dependencies lazily, so a stale module name
    survives until a real run reaches that line.
    """
    import ast
    import pathlib
    from importlib.util import resolve_name

    root = pathlib.Path(pymopsmap.__file__).parent
    for path in root.rglob("*.py"):
        package = ".".join(
            ("pymopsmap", *path.relative_to(root).with_suffix("").parts[:-1])
        )
        for node in ast.walk(ast.parse(path.read_text())):
            if not isinstance(node, ast.ImportFrom) or node.level == 0:
                continue
            name = "." * node.level + (node.module or "")
            importlib.import_module(resolve_name(name, package))
