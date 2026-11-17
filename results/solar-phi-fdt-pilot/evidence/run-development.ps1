$ErrorActionPreference = 'Stop'
$env:PYTHONPATH = 'D:\Open-Space\software\solar-phi-sunpy-20261117'
$python = 'C:\Users\Admin\source\repos\Open-Space\.venv\Scripts\python.exe'
$base = 'D:\Open-Space\campaigns\phi-fdt-4568cd0-metadata-terra61'
$out = 'D:\Open-Space\campaigns\phi-fdt-development-24ec48d-rootreview-20261117'
Set-Location 'D:\Open-Space\worktrees\solar-phi-fdt-20261117'
& $python results\solar-phi-fdt-pilot\development.py "$base\solo_L2_phi-fdt-blos_20260813T102509_V03_0648130504.fits" "$base\solo_L2_phi-fdt-blos_20260813T162509_V03_0648130506.fits" --holdout-meta-a "$base\solo_L2_phi-fdt-blos_20260815T162509_V03_0648150506.fits" --holdout-meta-b "$base\solo_L2_phi-fdt-blos_20260815T192509_V03_0648150507.fits" --run-id phi_fdt_dev_24ec48d_rootreview_20261117 --code-sha 24ec48d083e4214d5e2cb3eb5c17eaf3ccf169b2 --out "$out\development.json" *> "$out\development.log"
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
& $python results\solar-phi-fdt-pilot\prediction.py "$out\development.json" --out "$out\prediction.json" *> "$out\prediction.log"
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
