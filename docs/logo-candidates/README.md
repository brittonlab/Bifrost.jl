# Logo Candidates

Ten logo options for Bifrost.jl (issue #102). All candidates use the four Julia logo
colors — red `#CB3C33`, green `#389826`, purple `#9558B2`, blue `#4063D8` — and follow
Julia community conventions: flat geometric shapes, the three-dots motif where natural,
and a square 512×512 canvas suitable for `docs/src/assets/logo.svg`.

Each candidate is authored as hand-written SVG (the source of truth) with a 512 px PNG
render alongside it.

| File | Concept |
| --- | --- |
| `01-bridge-arc` | The Bifrost rainbow bridge as three nested arcs in Julia red, green, and purple, anchored by blue endpoints (Midgard and Asgard). |
| `02-julia-dots-fiber` | The classic Julia three-dots threaded by a blue optical fiber that arcs through them. |
| `03-birefringent-split` | Birefringence itself: one blue input ray splits inside a purple medium into ordinary (green) and extraordinary (red) rays. |
| `04-spun-helix` | Spun-fiber helix: three phase-shifted strands braiding along the fiber axis toward a blue end face. |
| `05-poincare-sphere` | The Poincaré sphere of polarization states with equator, meridian, and circular-polarization poles. |
| `06-fiber-cross-section` | Step-index fiber cross section: blue cladding, slightly elliptical purple core, and dashed fast/slow birefringence axes. |
| `07-gradient-b` | A bold letter B stroked with a continuous gradient through all four Julia colors. |
| `08-polarization-states` | Linear, circular, and elliptical polarization states arranged in the Julia three-dots layout. |
| `09-bent-fiber-strands` | A bent fiber path guiding three polarization strands — the package's core use case of propagation along curved geometry. |
| `10-dotted-bridge` | A rainbow bridge built entirely from Julia-colored dots. |

To regenerate the PNGs from the SVG sources:

```bash
for f in docs/logo-candidates/*.svg; do
    rsvg-convert -w 512 -h 512 "$f" -o "${f%.svg}.png"
done
```
