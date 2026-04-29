#!/usr/bin/env python3

# pyXLMS - pytest CONFIG
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--runslow", action="store_true", default=False, help="Run (all) slow tests."
    )
    parser.addoption(
        "--runext",
        action="store_true",
        default=False,
        help="Run (all) tests relying on external services.",
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow to run")
    config.addinivalue_line(
        "markers", "external: mark test as needing external service to run"
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow") and config.getoption("--runext"):
        # --runslow and --runext given in cli: do not skip tests
        return
    skip_slow = pytest.mark.skip(reason="needs --runslow option to run")
    skip_ext = pytest.mark.skip(reason="needs --runext option to run")
    for item in items:
        if "slow" in item.keywords:
            if not config.getoption("--runslow"):
                item.add_marker(skip_slow)
        elif "external" in item.keywords:
            if not config.getoption("--runext"):
                item.add_marker(skip_ext)
