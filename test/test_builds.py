# SPDX-FileCopyrightText: 2026 GFZ Helmholtz Centre for Geosciences
#
# SPDX-License-Identifier: BSD-3-Clause

# /// script
# dependencies = [
#   "rich",
# ]
# ///

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

from rich.console import Console

console = Console()

parser = argparse.ArgumentParser(
    description="Run the VERB4D build matrix tests."
)
parser.add_argument(
    "-l",
    "--openblas-libpath",
    default=None,
    help=(
        "Path to libopenblas.so to use for OpenBLAS-linked targets. "
        "If omitted, GetBlasType.cmake's default $HOME-based lookup is used."
    ),
)
args = parser.parse_args()

_openblas_libpath_override = (
    {"OPENBLAS_LIBPATH": args.openblas_libpath}
    if args.openblas_libpath
    else {}
)

default_target = {}  # use defaults

open_blas_target = {"BLAS_TYPE": "OpenBLAS"} | _openblas_libpath_override

data_assimilation_target = {
    "BLAS_TYPE": "OpenBLAS",
    "DATA_ASSIMILATION": "ON",
} | _openblas_libpath_override

data_assimilation_debug_target = {
    "BLAS_TYPE": "OpenBLAS",
    "DATA_ASSIMILATION": "ON",
    "DATA_ASSIMILATION_DEBUG": "ON",
} | _openblas_libpath_override

ppfv_target = {"USE_PPFV": "ON"}

python_bindings_target = {"PYTHON_BINDINGS": "ON"}

ppfv_python_bindings_target = {
    "USE_PPFV": "ON",
    "PYTHON_BINDINGS": "ON",
}

fast_convection_target = {"FAST_CONVECTION": "ON"}

save_psd_lost_target = {"SAVE_PSD_LOST_CONV": "ON"}

lu_caching_target = {"LU_CACHING": "ON"}

all_enabled_target = (
    default_target
    | open_blas_target
    | data_assimilation_target
    | data_assimilation_debug_target
    | ppfv_target
    | python_bindings_target
    | ppfv_python_bindings_target
    | fast_convection_target
    | save_psd_lost_target
    | lu_caching_target
)


test_targets = [
    ("default_target", default_target),
    ("open_blas_target", open_blas_target),
    ("data_assimilation_target", data_assimilation_target),
    ("data_assimilation_debug_target", data_assimilation_debug_target),
    ("ppfv_target", ppfv_target),
    ("python_bindings_target", python_bindings_target),
    ("ppfv_python_bindings_target", ppfv_python_bindings_target),
    ("fast_convection_target", fast_convection_target),
    ("save_psd_lost_target", save_psd_lost_target),
    ("lu_caching_target", lu_caching_target),
    ("all_enabled_target", all_enabled_target),
]

# uncomment to limit tests
# test_targets = [ppfv_python_bindings_target]

def run_build(config: dict, config_name: str, base_path: Path):
    build_dir = base_path / f"test_build_{config_name}"
    log_dir = base_path / "logs"

    if build_dir.exists():
        shutil.rmtree(build_dir)
    build_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(exist_ok=True)

    cmake_flags = [f"-D{k}={v}" for k, v in config.items()]
    cmake_flags.append("-DCMAKE_BUILD_TYPE=Release")

    # Run Configure
    conf = subprocess.run(  # noqa: PLW1510
        ["cmake", str(base_path.parent), "-B", str(build_dir)] + cmake_flags,
        capture_output=True,
        text=True,
    )
    if conf.returncode != 0:
        (log_dir / f"config_fail_{config_name}.log").write_text(conf.stderr)
        return False

    # Run Build
    build = subprocess.run(  # noqa: PLW1510
        ["cmake", "--build", str(build_dir), "--parallel", "4"],
        capture_output=True,
        text=True,
    )
    if build.returncode != 0:
        (log_dir / f"build_fail_{config_name}.log").write_text(build.stderr)
        return False

    return True


root_path = Path(__file__).parent.resolve()

console.print(
    "[bold red]Cleaning up old build directories and logs...[/bold red]"
)
# Find all directories matching test_build_* and remove them
for old_dir in root_path.glob("test_build_*"):
    if old_dir.is_dir():
        shutil.rmtree(old_dir)

# Remove old logs
old_logs = root_path / "logs"
if old_logs.exists():
    shutil.rmtree(old_logs)

console.print(
    f"[bold blue]Starting Matrix Build for {len(test_targets)} configurations...[/bold blue]\n"
)

failed_configs = []

with console.status("", spinner="dots2") as status:
    for i, (config_name, config) in enumerate(test_targets):
        # Format the current options for the display
        if len(config) == 0:
            opt_str = "default"
        else:
            opt_str = ", ".join(
                [f"[cyan]{k}[/]={v}" for k, v in config.items()]
            )
        status.update(
            f"[bold yellow]Config {config_name}:[/bold yellow] {opt_str}"
        )

        success = run_build(config, config_name, root_path)
        if not success:
            failed_configs.append(config_name)

        # Log result above the sticky status bar
        status_icon = "[green]PASS[/green]" if success else "[red]FAIL[/red]"
        console.log(f"{status_icon} Config {config_name}: {opt_str}")

if failed_configs:
    console.print(
        f"\n[bold red]{len(failed_configs)}/{len(test_targets)} configuration(s) failed:[/bold red] "
        + ", ".join(failed_configs)
    )
    sys.exit(1)

console.print(
    f"\n[bold green]All {len(test_targets)} configurations passed.[/bold green]"
)
