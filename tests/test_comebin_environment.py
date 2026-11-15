import os
import subprocess
import textwrap
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
COMEBIN_ENV = ROOT / "workflow" / "envs" / "comebin.yaml"
BINNING_RULES = ROOT / "workflow" / "rules" / "binning.smk"


def test_comebin_environment_pins_verified_cuda_build():
    assert yaml.safe_load(COMEBIN_ENV.read_text(encoding="utf-8")) == {
        "channels": ["pytorch", "conda-forge", "bioconda", "nodefaults"],
        "dependencies": [
            "comebin=1.0.4",
            "pytorch=1.10.2=py3.7_cuda11.3_cudnn8.2.0_0",
            "cudatoolkit=11.3",
            "mkl=2022",
            "protobuf=3.19.4",
        ],
    }


def test_comebin_rule_keeps_environment_reference_and_native_device_selection():
    comebin_rule = BINNING_RULES.read_text(encoding="utf-8").split(
        "\n\nrule metabinner_depth:", 1
    )[0]
    assert 'conda:\n        "../envs/comebin.yaml"' in comebin_rule
    assert "CUDA_VISIBLE_DEVICES" not in comebin_rule
    assert "--device" not in comebin_rule


def _supplied_environment_python():
    prefix = os.environ.get("MEEW_COMEBIN_ENV_PREFIX")
    if not prefix:
        pytest.skip("set MEEW_COMEBIN_ENV_PREFIX to test a solved COMEBin environment")
    python = Path(prefix) / "bin" / "python"
    if not python.is_file():
        pytest.fail(f"MEEW_COMEBIN_ENV_PREFIX has no executable Python: {python}")
    return python


def _run_smoke(python, script, tmp_path, visible_devices):
    env = os.environ.copy()
    env.update(
        {
            "CUDA_VISIBLE_DEVICES": visible_devices,
            "PYTHONDONTWRITEBYTECODE": "1",
            "TMPDIR": str(tmp_path),
            "XDG_CACHE_HOME": str(tmp_path / "cache"),
            "MPLCONFIGDIR": str(tmp_path / "matplotlib"),
        }
    )
    return subprocess.run(
        [str(python), "-c", textwrap.dedent(script)],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )


def test_supplied_comebin_environment_retains_cpu_fallback(tmp_path):
    python = _supplied_environment_python()
    result = _run_smoke(
        python,
        """
        import inspect
        import os
        import sys

        sys.path.insert(0, os.path.join(sys.prefix, "bin", "COMEBin"))
        import torch
        import train_CLmodel
        from models.mlp import EmbeddingNet

        source = inspect.getsource(train_CLmodel.train_CLmodel)
        assert "if torch.cuda.is_available():" in source
        assert "torch.device('cuda')" in source
        assert "torch.device('cpu')" in source
        assert not torch.cuda.is_available()
        assert torch.cuda.device_count() == 0

        model = EmbeddingNet(8, 4, [16, 16], [0.0])
        optimizer = torch.optim.AdamW(model.parameters(), lr=0.001)
        output = model(torch.randn(4, 8))
        loss = output.square().mean()
        optimizer.zero_grad()
        loss.backward()
        assert torch.isfinite(output).all()
        assert all(
            parameter.grad is None or torch.isfinite(parameter.grad).all()
            for parameter in model.parameters()
        )
        optimizer.step()
        print("COMEBin CPU fallback smoke passed")
        """,
        tmp_path,
        "",
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_supplied_comebin_environment_cuda_smoke(tmp_path):
    python = _supplied_environment_python()
    device = os.environ.get("MEEW_COMEBIN_CUDA_DEVICE")
    if device is None:
        pytest.skip("set MEEW_COMEBIN_CUDA_DEVICE to a suitable free GPU index")
    result = _run_smoke(
        python,
        """
        import os
        import sys

        sys.path.insert(0, os.path.join(sys.prefix, "bin", "COMEBin"))
        import torch
        import torch.nn as nn
        from torch.cuda.amp import GradScaler, autocast
        import simclr
        import train_CLmodel
        from models.mlp import EmbeddingNet
        from models.mlp2 import EmbeddingNet as CoverageEmbeddingNet

        assert torch.__version__ == "1.10.2"
        assert torch.version.cuda == "11.3"
        assert torch.cuda.is_available()
        assert torch.cuda.device_count() == 1
        device = torch.device("cuda:0")

        coverage_model = EmbeddingNet(4, 3, [8, 8], [0.0], actn=nn.LeakyReLU())
        model = CoverageEmbeddingNet(
            139,
            6,
            [16, 16],
            [0.0],
            actn=nn.LeakyReLU(),
            cov_model=coverage_model,
        ).to(device)
        optimizer = torch.optim.AdamW(model.parameters(), lr=0.001)
        scaler = GradScaler(enabled=True, init_scale=1.0)
        optimizer.zero_grad()
        with autocast(enabled=True):
            output, coverage_output = model(
                torch.randn(8, 136, device=device),
                torch.randn(8, 4, device=device),
            )
            loss = output.float().square().mean()
            loss = loss + coverage_output.float().square().mean()
        scaler.scale(loss).backward()
        scaler.unscale_(optimizer)
        gradients = [p.grad for p in model.parameters() if p.grad is not None]
        assert gradients
        assert torch.isfinite(output).all()
        assert torch.isfinite(coverage_output).all()
        assert torch.isfinite(loss)
        assert all(torch.isfinite(gradient).all() for gradient in gradients)
        scaler.step(optimizer)
        scaler.update()
        torch.cuda.synchronize()
        print("COMEBin CUDA model smoke passed", torch.cuda.get_device_name(0))
        """,
        tmp_path,
        device,
    )
    assert result.returncode == 0, result.stdout + result.stderr
