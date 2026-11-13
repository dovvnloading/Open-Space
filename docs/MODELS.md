# Model-assisted research

The analyst receives a numerical evidence packet from one spectrum and can request up to five targeted measurements. The worker computes those measurements against the original input rows. The reviewer then receives the evidence, the analyst's explanation and the measured outcomes. A reviewer can identify a competing hypothesis, describe an observation that would discriminate between the two, or conclude that the evidence is insufficient. Everything a model writes is recorded as a proposal or judgment attributed to that model. The local measurements retain their own values and provenance.

Use the default deterministic mode to run the archive, detector and validation controls offline. For model mode, run a chat-completions server on your compute or choose a hosted endpoint that supports JSON object responses and the `max_completion_tokens` parameter. A self-hosted vLLM deployment requires a served model with a chat template. Select the model identifier exposed by your server. You can assign independent models to the two roles. [vLLM server documentation](https://docs.vllm.ai/en/latest/serving/openai_compatible_server/)

```sh
export OPEN_SPACE_MODEL_BASE_URL=http://127.0.0.1:8000/v1
export OPEN_SPACE_ANALYST_MODEL=your-served-model-id
export OPEN_SPACE_REVIEWER_MODEL=your-reviewer-model-id
export OPEN_SPACE_MAX_TOKENS=1600

python -m open_space enqueue --models data/sdss-6497-56329-0139.json
python -m open_space work --max-jobs 1 --max-calls 2
```

| Setting | Behavior |
| --- | --- |
| `OPEN_SPACE_MODEL_BASE_URL` | Required in model mode. The client appends `/chat/completions`. Use HTTPS for hosted endpoints. |
| `OPEN_SPACE_ANALYST_MODEL` | Required model identifier. |
| `OPEN_SPACE_REVIEWER_MODEL` | Defaults to the analyst model. Differences between replies from one or more models do not demonstrate independent observational replication. |
| `OPEN_SPACE_MAX_TOKENS` | Requested maximum completion tokens per call, including any reasoning tokens the server counts. Defaults to 1,600. Accepts integers from 128 to 8,192. Input tokens and endpoint pricing are separate. |
| `OPEN_SPACE_API_KEY` | Optional bearer key for the selected endpoint. Supply it through your environment or a secret manager. The client reads it when making the request. |

## What each model sees

The client assigns IDs in the original detector order. `F0001` is the first feature in the saved analysis. It sends at most 24 features, prioritizing those outside the bundled catalog, then ordering each group by descending score with original index as a tiebreaker. The packet records total and omitted counts. Each feature includes the signed contrast, direction, both wavelength frames, local pixel spacing, matching tolerance and up to four catalog matches. The quality summary and counts cover the full measurement. The client sends the packet and its role instructions to the configured endpoint. The worker includes the validated observation ID and source kind. If another caller supplies no source context, the client marks the origin unspecified.

| Agent | Required output | Available decision |
| --- | --- | --- |
| Analyst | Summary, hypotheses with their predictions and falsifiers, requested measurements | Request `mask_neighborhood` or `continuum_sensitivity` for a supplied feature ID. An empty request list is valid. |
| Reviewer | Summary and feature-specific assessments, each with a challenge and next test | Record `explained`, `inconclusive` or `follow_up`. |

Both response schemas limit each list to five entries. Validation rejects unknown IDs, duplicate measurement requests, unknown fields, unsupported checks, incomplete responses and excessive text. Model content is parsed as JSON. Only the two declared check names can select code paths in the worker. When serving a hosted model, configure its retention and data settings to suit the observations you intend to provide.

## Reproducing a run

Each returned role contains a `review` object with the protocol version, role, endpoint, requested and returned model IDs, exact system and user messages, completion cap, provider response ID, reported token counts when available, and SHA-256 digests of the request, response and validated output. Hashes use UTF-8 JSON with sorted keys, no extra whitespace and unescaped Unicode. Preserve the weights revision, chat template, quantization, inference-server version and decoding configuration separately for self-hosted experiments. A recorded prompt permits an audit even when a provider changes weights or sampling behavior.

A role makes one POST with a 60-second socket timeout and no automatic retry or redirect. The client accepts at most 128,000 response bytes and a 64,000-byte evidence packet. The worker counts model calls across attempts within each work invocation. The `--max-calls` budget applies across all jobs in that invocation. A new invocation has its own budget. The default database is `runs/ledger.sqlite`. Use `python -m open_space --db PATH ...` to select a different ledger. A lost connection after the server accepted a request may still consume tokens. Check the server's request logs before rerunning a failed role. The test suite uses a local HTTP server to exercise both roles without sending observations to an external service.
