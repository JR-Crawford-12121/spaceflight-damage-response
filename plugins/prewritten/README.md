# Prewritten plugins (future AutoBioScout)

This directory holds **analysis plugins** that a future **AutoBioScout** agent may invoke after the base pipeline runs.

Planned loop (not implemented here):

1. Read **`outputs/logs/initial_analysis_summary.json`** produced by `main.py`.
2. Identify uncertainty (sparse overlaps, conflicting pathways, missing mappings).
3. Select a follow-up analysis implemented as a **prewritten plugin** (e.g., composition-aware check, robustness sweep).
4. Execute the plugin in a controlled subprocess or import path.
5. Append results and decisions to an **audit log** under `plugins/generated/` or `outputs/logs/`.

Keep plugins deterministic and **free of chat/completions requirements** unless explicitly configured later.
