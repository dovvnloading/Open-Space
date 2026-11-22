"""Restrict archive access and survey computation to the authorized Actions workflow."""
import os


def require_survey_workflow():
    event = os.environ.get("GITHUB_EVENT_NAME")
    branch = "refs/heads/survey/harps-rvbank-companions"
    workflow = os.environ.get("GITHUB_WORKFLOW_REF", "").partition("@")[0]
    expected = os.environ.get("GITHUB_REPOSITORY", "") + "/.github/workflows/harps-rvbank-survey.yml"
    permitted = event == "workflow_dispatch" or (event == "push" and os.environ.get("GITHUB_REF") == branch)
    if os.environ.get("GITHUB_ACTIONS") != "true" or workflow != expected or not permitted:
        raise SystemExit("Run the survey in the HARPS Actions workflow.")
