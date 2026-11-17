$ErrorActionPreference = 'Stop'
$env:PYTHONPATH='D:\Open-Space\software\solar-phi-sunpy-20261117'
$python='C:\Users\Admin\source\repos\Open-Space\.venv\Scripts\python.exe'
$base='D:\Open-Space\campaigns\phi-fdt-4568cd0-metadata-terra61'
$out='D:\Open-Space\campaigns\phi-fdt-holdout-24ec48d-973bd0f-20261117'
$prediction='D:\Open-Space\campaigns\phi-fdt-development-24ec48d-rootreview-20261117\prediction.json'
Set-Location 'D:\Open-Space\worktrees\solar-phi-fdt-20261117'
& $python results\solar-phi-fdt-pilot\validate_repeat.py --prediction $prediction --authorization "$out\holdout-authorization.json" --holdout-a "$base\solo_L2_phi-fdt-blos_20260815T162509_V03_0648150506.fits" --holdout-b "$base\solo_L2_phi-fdt-blos_20260815T192509_V03_0648150507.fits" --out "$out\primary-result.json" *> "$out\primary.log"
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
