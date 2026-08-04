from __future__ import annotations

from figure_config import resolve_external_sources


def test_cli_values_override_environment(monkeypatch, tmp_path):
    env_root = tmp_path / "environment"
    cli_root = tmp_path / "cli"
    monkeypatch.setenv("STRIDE_OUTPUTS_ROOT", str(env_root))
    resolved = resolve_external_sources(outputs_root=cli_root)
    assert resolved.outputs_root == cli_root.resolve()


def test_environment_values_override_relative_defaults(monkeypatch, tmp_path):
    values = {
        "STRIDE_PRODUCTION_ROOT": tmp_path / "production",
        "STRIDE_DOWNSTREAM_ROOT": tmp_path / "downstream",
        "STRIDE_OUTPUTS_ROOT": tmp_path / "outputs",
        "STRIDE_MANUSCRIPT_PATH": tmp_path / "manuscript.md",
        "STRIDE_V9A_ARTIFACTS": tmp_path / "v9a",
        "STRIDE_V9B_ARTIFACTS": tmp_path / "v9b",
    }
    for name, value in values.items():
        monkeypatch.setenv(name, str(value))
    resolved = resolve_external_sources()
    assert resolved.production_root == values["STRIDE_PRODUCTION_ROOT"].resolve()
    assert resolved.downstream_root == values["STRIDE_DOWNSTREAM_ROOT"].resolve()
    assert resolved.outputs_root == values["STRIDE_OUTPUTS_ROOT"].resolve()
    assert resolved.manuscript_path == values["STRIDE_MANUSCRIPT_PATH"].resolve()
    assert resolved.v9a_artifacts == values["STRIDE_V9A_ARTIFACTS"].resolve()
    assert resolved.v9b_artifacts == values["STRIDE_V9B_ARTIFACTS"].resolve()


def test_v9a_and_v9b_defaults_derive_from_repository_roots(tmp_path):
    production = tmp_path / "production"
    downstream = tmp_path / "downstream"
    resolved = resolve_external_sources(
        production_root=production, downstream_root=downstream
    )
    assert resolved.v9a_artifacts == (
        production / "validation" / "artifacts" / "v9a"
    ).resolve()
    assert resolved.v9b_artifacts == (downstream / "validation_artifacts").resolve()
