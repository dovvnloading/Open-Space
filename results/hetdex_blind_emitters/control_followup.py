"""Apply the same source vetting to reconstructed matched and raw-shuffle spectra."""
from nuisance_control import records_for_vet, raw_subsample_for_vet, run_control
from vetting import vet


def control_followup(catalog, output):
    run_control(catalog, output)
    vet(records_for_vet(catalog, output), catalog, output / 'matched-shuffle-vet')
    vet(raw_subsample_for_vet(catalog, output, limit=1024), catalog, output / 'shuffle-vet')
