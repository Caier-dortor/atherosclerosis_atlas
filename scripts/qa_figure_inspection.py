"""
SCI Figure QA Inspector — checks all 7 figures against publication standards.
"""
from pathlib import Path
import re
import sys
# Fix Windows console encoding
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

RES_DIR = Path("D:/openclaw_workspace/atherosclerosis_atlas/results")
FIG_DIRS = ['fig1', 'fig2', 'fig3', 'fig4', 'fig5', 'fig6', 'fig7']

print("=" * 70)
print("SCI FIGURE QA INSPECTION REPORT")
print("=" * 70)

all_issues = []
all_passes = []
all_warnings = []

for fig_dir in FIG_DIRS:
    fig_path = RES_DIR / fig_dir
    png_files = list(fig_path.glob("*.png"))
    svg_files = list(fig_path.glob("*.svg"))

    if not png_files:
        all_issues.append(f"{fig_dir}: No PNG found")
        continue

    png = png_files[0]
    svg = svg_files[0] if svg_files else None

    print(f"\n{'='*70}")
    print(f"  {fig_dir.upper()}: {png.name}")
    print(f"{'='*70}")

    # === 1. PNG DPI ===
    try:
        from PIL import Image
        img = Image.open(png)
        dpi = img.info.get('dpi', (None, None))
        w_px, h_px = img.size

        if dpi[0] and dpi[0] >= 290:  # allow 10 DPI tolerance
            w_in = w_px / dpi[0]
            all_passes.append(f"{fig_dir}: DPI={dpi[0]:.0f} (OK)")
            print(f"  PASS  DPI: {dpi[0]:.0f} ({w_px}x{h_px} px, {w_in:.1f} in wide)")
        elif dpi[0]:
            all_issues.append(f"{fig_dir}: DPI={dpi[0]:.0f} < 300")
            print(f"  FAIL  DPI={dpi[0]:.0f} < 300")
        else:
            all_warnings.append(f"{fig_dir}: No DPI metadata")
            print(f"  WARN  No DPI metadata ({w_px}x{h_px} px)")
    except ImportError:
        print("  SKIP  PIL not available")

    # === 2. SVG analysis ===
    if svg and svg.exists():
        svg_text = svg.read_text(encoding='utf-8', errors='replace')
        svg_size = svg.stat().st_size / 1024
        png_size = png.stat().st_size / 1024

        # Check for real <text> elements (not path-based text)
        text_tags = re.findall(r'<text[ >]', svg_text)
        has_text_tags = len(text_tags) > 0

        # Check for path-based text (glyphs)
        glyph_paths = re.findall(r'id="([^"]*PSMT[^"]*)"', svg_text)
        has_glyph_paths = len(glyph_paths) > 0

        # Check font-family in styles or text
        font_families = set(re.findall(r'font-family:([^;"]+)[;"]', svg_text))

        if has_text_tags:
            all_passes.append(f"{fig_dir}: {len(text_tags)} editable <text> elements")
            print(f"  PASS  Editable text: {len(text_tags)} <text> elements")
        elif has_glyph_paths:
            all_issues.append(f"{fig_dir}: Text is path-based (un-editable), {len(glyph_paths)} glyph paths")
            print(f"  FAIL  Text is paths, NOT editable ({len(glyph_paths)} glyph paths)")
        else:
            all_warnings.append(f"{fig_dir}: No text elements found")
            print(f"  WARN  No text elements detected")

        # Font check
        if font_families:
            has_tnr = any('Times' in f for f in font_families)
            if has_tnr:
                all_passes.append(f"{fig_dir}: Times font family detected")
                print(f"  PASS  Font: Times-family detected")
            else:
                print(f"  INFO  Font families: {font_families}")
        else:
            print(f"  INFO  Font: embedded as paths (Times New Roman glyphs)")

        # === 3. Color analysis ===
        hex_colors = set(re.findall(r'#[0-9A-Fa-f]{6}', svg_text))

        # CB palette check (case-insensitive)
        cb_map = {'#D55E00': 'carotid', '#0072B2': 'coronary', '#009E73': 'femoral'}
        cb_found = [c for c in hex_colors if c.upper() in [k.upper() for k in cb_map]]
        if len(cb_found) >= 2:
            all_passes.append(f"{fig_dir}: CB palette detected ({len(cb_found)}/3)")
            print(f"  PASS  CB palette: {len(cb_found)}/3 detected")
        elif len(cb_found) == 1:
            all_warnings.append(f"{fig_dir}: Only 1/3 CB colors")
            print(f"  WARN  CB palette: only {len(cb_found)}/3")

        # Check for RdBu_r colormap (heatmap safety)
        if 'RdBu' in svg_text or 'rdbu' in svg_text.lower():
            print(f"  INFO  Heatmap colormap: RdBu_r (colorblind-safe)")
        if 'YlOrRd' in svg_text:
            print(f"  INFO  Heatmap colormap: YlOrRd")

        # === 4. Panel structure ===
        panel_letters = set(re.findall(r'>([A-F])</text>', svg_text))
        if not panel_letters:
            # Try to find in fig.text positioning (rasterized)
            panel_letters = set(re.findall(r'\(([A-F])\)', svg_text))

        if panel_letters:
            print(f"  PASS  Panels: {sorted(panel_letters)}")
        elif has_glyph_paths:
            print(f"  INFO  Panels likely embedded in path text")

        # === 5. Statistical annotations ===
        sig_markers = ['***', '**', '*  ', 'ns  ', 'p=', 'rs=', 'H=', 'KW']
        found = [s for s in sig_markers if s in svg_text]
        if found:
            print(f"  PASS  Stats markers in SVG: {found}")

        # === 6. Figure dimensions ===
        if dpi[0] and dpi[0] > 0:
            w_in = w_px / dpi[0]
            h_in = h_px / dpi[0]
            print(f"  INFO  Size: {w_in:.1f}x{h_in:.1f} in, PNG={png_size:.0f}KB, SVG={svg_size:.0f}KB")
            if w_in > 15:
                all_warnings.append(f"{fig_dir}: Very wide ({w_in:.0f}in) — check journal max width")

        # === 7. Check fig.text() vs axes titles (panel labels) ===
        # Look for fig.text positioned labels (bold panel letters)
        figtext_matches = re.findall(r'font-weight:bold[^>]*>([A-F])<', svg_text)
        if figtext_matches:
            print(f"  INFO  Panel labels via fig.text(): {figtext_matches}")

    else:
        all_issues.append(f"{fig_dir}: No SVG")
        print(f"  FAIL  No SVG found")

# ============================================================
# Cross-figure consistency
# ============================================================
print(f"\n{'='*70}")
print("CROSS-FIGURE CONSISTENCY")
print(f"{'='*70}")

# Check all SVG sizes
svg_sizes = {}
for fig_dir in FIG_DIRS:
    svg_files = list((RES_DIR / fig_dir).glob("*.svg"))
    if svg_files:
        svg_sizes[fig_dir] = svg_files[0].stat().st_size / 1024

if svg_sizes:
    avg_svg = sum(svg_sizes.values()) / len(svg_sizes)
    for fig, sz in svg_sizes.items():
        if sz < avg_svg * 0.3:
            all_warnings.append(f"{fig}: SVG unusually small ({sz:.0f}KB vs avg {avg_svg:.0f}KB)")

print(f"SVG sizes: {', '.join(f'{k}={v:.0f}KB' for k, v in sorted(svg_sizes.items()))}")
print(f"Average SVG: {avg_svg:.0f}KB")

# ============================================================
# Summary
# ============================================================
print(f"\n{'='*70}")
print("QA SUMMARY")
print(f"{'='*70}")

print(f"\nPASSES ({len(all_passes)}):")
for p in all_passes:
    print(f"  + {p}")

print(f"\nISSUES ({len(all_issues)}):")
for i in all_issues:
    print(f"  ! {i}")

if all_warnings:
    print(f"\nWARNINGS ({len(all_warnings)}):")
    for w in all_warnings:
        print(f"  ~ {w}")

# Critical findings
critical = [i for i in all_issues if 'text is path' in i.lower() or 'un-editable' in i.lower()]
if critical:
    print(f"\n>> CRITICAL: Text editability ({len(critical)} figs)")
    print(f"   Root cause: 'svg.fonttype': 'none' is deprecated in matplotlib 3.10")
    print(f"   Fix: Change to 'svg.fonttype': 'svgfont' in all scripts")

print(f"\nScore: {len(all_passes)}/{len(all_passes)+len(all_issues)} checks passed")
print(f"{'='*70}")
