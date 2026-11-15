// Technical note. Build with `make note` (writes docs/technical-note.pdf).
// Every number comes from reports/results.json, written by `make reproduce`.

#let r = json("/reports/results.json")
#let n(x, d: 2) = str(calc.round(x, digits: d))

#set document(title: "Density and E-field coupling in the dayside aurora", author: "Renato Valentim", date: datetime(year: 2026, month: 11, day: 15))
#set page(paper: "us-letter", margin: (x: 0.8in, y: 0.7in), numbering: "1 / 1")
#set text(font: "Libertinus Serif", size: 10pt)
#set par(justify: true, leading: 0.55em)
#show heading.where(level: 1): it => block(above: 1.0em, below: 0.5em, text(size: 10.5pt, weight: "bold", it.body))
#show figure.caption: set text(size: 8.5pt)

#align(center)[
  #text(size: 14pt, weight: "bold")[Density and electric-field coupling in the dayside aurora: \ a reproducible re-analysis of TRICE-2 High data]
  #v(0.3em)
  Renato Valentim \
  #text(size: 9pt)[Research supervised by Prof. J. LaBelle, Dartmouth College (2021-2022). Technical note, November 2026. \ Code and data: #link("https://github.com/valentim7/plasma_physics")]
]

#v(0.4em)
#block(inset: (x: 1.2em))[
  #text(size: 9.5pt)[
    *Summary.* We test whether electron-density fluctuations (from the traced plasma frequency) and electric-field fluctuations measured on the TRICE-2 High rocket follow the Boltzmann relation. The 2022 MATLAB analysis is reproduced exactly, the sun spikes are now removed automatically, and the coupling is tested directly with cross-spectral coherence. Coherence at 10-40 Hz is small but significant over the flight (permutation test, $p = #n(r.coh_p, d: 3)$) and is concentrated in #r.n_sig_blocks of #r.n_blocks blocks of about 28 s, all between #calc.round(r.sig_span.at(0)) and #calc.round(r.sig_span.at(1)) s, where wave activity is strongest. The main limitation is the density itself: before about 500 s it is at its noise floor in this band. An automatic plasma-frequency trace from the raw HF spectra is the natural next step.
  ]
]

= 1. Question and data

The Boltzmann relation links electron-density fluctuations to the electric potential, $delta n slash n = e delta phi slash k T_e$. If it holds, density and electric-field fluctuations must vary together. The density comes from the plasma frequency traced on HF spectrograms ($n_e = (f_"pe" slash 8.98)^2$, $f_"pe"$ in kHz, $n_e$ in cm#super[-3]); the trace was found semi-automatically and then corrected point by point by hand. The electric field is a single channel sampled at 2500 Hz, in what appear to be raw instrument counts. It shows a large oscillation with a period of about 1.84 s and a short "sun spike" every half of that period, alternating in sign. This is consistent with a spinning double-probe instrument; the geometry is not documented in our files and is inferred from the signal.

= 2. Methods

*Reproduction.* The MATLAB slope analysis was ported to Python, quirks included. All #r.reproduced windows match the original output (maximum difference $#calc.round(r.reproduction_max_diff * 1e7) times 10^(-7)$, the rounding of the original text files).

*Sun spikes.* Spikes are detected from sample-to-sample jumps far above a robust (median absolute deviation) threshold: #r.spikes spikes, median spacing #n(r.spike_spacing, d: 3) s. Only gaps of 0.75-0.9 s between consecutive spikes are kept, so a missed spike cannot leak into a segment. This gives #r.segments clean segments of 0.75 s, against 64 windows picked by hand in 2022.

*Coherence.* The E field is low-pass filtered (zero phase, 60 Hz) and interpolated onto the density time grid; both series are detrended and Hann-windowed per segment, and the magnitude-squared coherence is averaged over segments. Because the spikes alternate in sign, which suggests the sensor flips orientation every half spin, each E segment is multiplied by the polarity of the preceding spike; without this, real correlation cancels (coherence #n(r.coh_raw, d: 4) without and #n(r.coh_corrected, d: 4) with the correction). Significance comes from #r.n_perm permutations that re-pair density and E segments at random, over the whole flight and per block of 30 segments (Bonferroni-corrected).

#figure(
  image("/reports/figures/coherence_blocks.png", width: 92%),
  caption: [Coherence between density and E field at 10-40 Hz per block of 30 clean segments. Gray: 95th percentile of the coherence obtained by shuffling the segment pairings.],
)

= 3. Results

- *Coupling exists but is localized.* Over the flight the coherence is #n(r.coh_corrected, d: 4) against a null mean of #n(r.coh_null_mean, d: 4) ($p = #n(r.coh_p, d: 3)$). It is significant in #r.n_sig_blocks of #r.n_blocks blocks, all between about #calc.round(r.sig_span.at(0)) and #calc.round(r.sig_span.at(1)) s, with block values up to #n(r.max_block_coh). This matches the 2022 impression that the relation appears only in some intervals; even there it is weak (about 10% shared variance).

- *Coherence follows wave activity (exploratory).* Per block, coherence rises with E-field power at 10-40 Hz (Spearman $rho = #n(r.rho_e)$, $p = #n(r.p_e, d: 3)$) and with density signal-to-noise ($rho = #n(r.rho_n)$). The strongest blocks appear to be where the rocket is inside the aurora. Two readings fit: the coupling is physically stronger there, or it is only detectable when both signals are strong. Before about 490 s the density power at 10-40 Hz equals its 60-75 Hz floor, so the data cannot separate the two.

- *Spectral slopes are a fragile test.* In 500-815 s, E-field slope minus density slope is #n(r.auto_diff_mean), 95% CI [#n(r.auto_diff_ci.at(0)), #n(r.auto_diff_ci.at(1))] (moving-block bootstrap). On the 2022 hand-picked windows it ranges from #n(r.band_diff_min) to #n(r.band_diff_max) depending on the fit band. A working expectation of +2 follows from $E = -nabla phi$ and frozen-in structures, but it assumes constant $T_e$ and a measured component along the wave vector, which a single rotating channel would not guarantee.

#figure(
  image("/reports/figures/coherence_vs_power.png", width: 62%),
  caption: [Block coherence against E-field power at 10-40 Hz (relative to the quietest block). Filled: significant after Bonferroni correction.],
)

= 4. Next step: an automatic plasma-frequency trace

The density noise floor sets the limit of every result above, so the next step is to rebuild the trace from the raw HF spectra (about 50 million points for this flight, not in the public repository). Proposed method: (1) smooth each log spectrum and score every candidate frequency by the drop in power across it; (2) find the best path through time with the Viterbi algorithm, penalizing large jumps and staying below the upper-hybrid cutoff, which removes the need for manual corrections of isolated spikes and dips; (3) refine each point with a local step fit and attach an uncertainty from the forward-backward posterior. Validation: agreement with the existing hand-corrected trace (about 126,000 points) on a held-out part of the flight, and a lower density noise floor at 10-40 Hz. Coherence is computed only afterwards and is never used to tune the trace. The same pipeline would then apply to TRICE-2 Low, CAPER-2 and other rocket data sets.

= 5. Limitations

The instrument geometry (spin, probe configuration, calibration, frequency response) is inferred from the signal and needs confirmation; with a single channel the field cannot be despun, and no instrument response is applied. The spike detector is conservative (some small jumps are flagged too), which discards segments but never lets a spike in. Coherence blocks have a fixed number of segments, so their duration varies (28 to 109 s). The block-level analysis of Section 3 was not planned in advance and should be read as exploratory. Two of the original MATLAB scripts carry a "TRICE2-LOW" header; the data is believed to be from TRICE-2 High, to be confirmed.

#v(0.3em)
#text(size: 8.5pt)[*Reproducibility.* `make all` sets up Python 3.12, verifies data checksums, runs 11 tests (including the exact MATLAB reproduction and synthetic checks of each method) and regenerates every number and figure in this note.]
