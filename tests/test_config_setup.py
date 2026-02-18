from pathlib import Path

from gollum import precomputed_spectrum as ps


def test_get_config_value_uses_template_fallback_without_writing_env(tmp_path, monkeypatch):
    config_env = tmp_path / "config.env"
    config_template = tmp_path / "config_template.env"
    config_template.write_text('PHOENIX = "~/template/path/"\n', encoding="utf-8")

    monkeypatch.setattr(ps, "CONFIG_ENV_PATH", config_env)
    monkeypatch.setattr(ps, "CONFIG_TEMPLATE_PATH", config_template)

    value = ps.get_config_value("PHOENIX")

    assert value == "~/template/path/"
    assert not config_env.exists()


def test_initialize_config_env_creates_file_from_template(tmp_path, monkeypatch):
    config_env = tmp_path / "config.env"
    config_template = tmp_path / "config_template.env"
    config_template.write_text('SonoraB = "~/template/bobcat/"\n', encoding="utf-8")

    monkeypatch.setattr(ps, "CONFIG_ENV_PATH", config_env)
    monkeypatch.setattr(ps, "CONFIG_TEMPLATE_PATH", config_template)

    out_path = ps.initialize_config_env()

    assert out_path == config_env
    assert config_env.exists()
    assert config_env.read_text(encoding="utf-8") == config_template.read_text(
        encoding="utf-8"
    )


def test_get_config_value_prefers_config_env_over_template(tmp_path, monkeypatch):
    config_env = tmp_path / "config.env"
    config_template = tmp_path / "config_template.env"
    config_env.write_text('coolTLUSTY = "~/env/path/"\n', encoding="utf-8")
    config_template.write_text('coolTLUSTY = "~/template/path/"\n', encoding="utf-8")

    monkeypatch.setattr(ps, "CONFIG_ENV_PATH", config_env)
    monkeypatch.setattr(ps, "CONFIG_TEMPLATE_PATH", config_template)

    assert ps.get_config_value("coolTLUSTY") == "~/env/path/"
