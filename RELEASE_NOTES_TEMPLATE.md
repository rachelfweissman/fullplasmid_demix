# PlasmidDemix v1.0.0 — desktop app

First standalone macOS release of the Plasmid Demix pipeline as a
double-clickable app.

## Install

1. Download `PlasmidDemix.dmg` below.
2. Open it and drag **PlasmidDemix** into your `Applications` folder.
3. **Important — one-time Terminal step:** macOS will say the app is
   *"damaged and can't be opened"* the first time. This is a false
   alarm caused by macOS quarantining unsigned downloads. Fix it with:
   ```bash
   xattr -cr /Applications/PlasmidDemix.app
   ```
   Then double-click the app normally. You only have to do this once.
4. First launch takes a few minutes while the bioinformatics tools
   (samtools, minimap2, flye, racon, medaka, …) install into
   `~/.plasmid_demix/`. After that, the app opens instantly.

## What's in this release

- Click-to-run macOS desktop app — no Python or conda setup needed.
- Live log panel showing pipeline progress.
- **Stop** button that cleanly kills any running tool (flye, medaka, etc.)
  if you need to abort partway through.
- Wraps the existing CLI (`fullPlasmidSeq_demix_RFW.py`); identical outputs.

## Why the macOS "damaged" warning?

The app is not signed with an Apple Developer ID ($99/year). On Apple
Silicon Macs running recent macOS versions, unsigned apps downloaded
from the internet are flagged with `com.apple.quarantine`, which now
shows up as the misleading "damaged" message. The `xattr -cr` command
above clears that flag.
