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
    *Summary.* We test whether electron-density fluctuations (from the traced plasma frequency) and electric-field fluctuations measured on the TRICE-2 High rocket follow the Boltzmann relation. The 2022 MATLAB analysis is reproduced exactly, the sun spikes are removed automatically, and the coupling is tested with the cross-spectrum. At 10-40 Hz the two signals are coupled in quadrature (90° phase shift, no in-phase part) in #r.n_sig_blocks of #r.n_blocks blocks of 30 clean segments, the phase expected when density follows the potential and the field is its gradient. The sign of the quadrature flips mid-flight. The coupling is found where wave activity is strong; elsewhere the density is at its noise floor, so an automatic plasma-frequency trace from the raw HF spectra is the natural next step.
  ]
]

= 1. Question and data

The Boltzmann relation links electron-density fluctuations to the electric potential, $delta n slash n = e delta phi slash k T_e$. If it holds, density and electric-field fluctuations must vary together. The density comes from the plasma frequency traced on HF spectrograms ($n_e = (f_"pe" slash 8.98)^2$, $f_"pe"$ in kHz, $n_e$ in cm#super[-3]); the trace was found semi-automatically and then corrected point by point by hand. The electric field is a single channel sampled at 2500 Hz, in what appear to be raw instrument counts. It shows a large oscillation with a period of about 1.84 s and a short "sun spike" every half of that period, alternating in sign. This is consistent with a spinning double-probe instrument; the geometry is not documented in our files and is inferred from the signal.

= 2. Methods

*Reproduction.* The MATLAB slope analysis was ported to Python, quirks included. All #r.reproduced windows match the original output (maximum difference $#calc.round(r.reproduction_max_diff * 1e7) times 10^(-7)$, the rounding of the original text files).

*Sun spikes.* Spikes are detected from sample-to-sample jumps far above a robust (median absolute deviation) threshold: #r.spikes spikes, median spacing #n(r.spike_spacing, d: 3) s. Only gaps of 0.75-0.9 s between consecutive spikes are kept, so a missed spike cannot leak into a segment. This gives #r.segments clean segments of 0.75 s, against 64 windows picked by hand in 2022.

*Coupling.* The E field is low-pass filtered (zero phase, 60 Hz) and interpolated onto the density time grid; both series are detrended and Hann-windowed per segment. The spikes alternate in sign, which suggests the sensor flips orientation every half spin, so each E segment is multiplied by the polarity of the preceding spike. For each segment we compute the band-integrated cross-spectrum divided by that segment's own power (a complex coherency: real part in phase, imaginary part in quadrature), and average it over blocks of 30 segments. Every segment therefore weighs the same. Significance comes from #r.n_perm shuffles of the segment pairings within each block, Bonferroni-corrected over #r.n_blocks blocks.

#figure(
  image("/reports/figures/quadrature_blocks.png", width: 90%),
  caption: [Mean coherency between density and E field at 10-40 Hz per block of 30 clean segments. Gray: 95% range obtained by shuffling the segment pairings.],
)

= 3. Results

- *Quadrature coupling in part of the flight.* #r.n_sig_blocks of #r.n_blocks blocks show significant quadrature and #r.n_sig_real show a significant in-phase part. The strongest block (#calc.round(r.strongest_span.at(0))-#calc.round(r.strongest_span.at(1)) s) has a quadrature of #n(r.strongest_quad), #n(calc.abs(r.strongest_z), d: 1) standard deviations from shuffled pairings. An interference common to both channels would appear in the in-phase part, which is not seen.

- *The sign flips.* The quadrature is positive in #calc.round(r.pos_span.at(0))-#calc.round(r.pos_span.at(1)) s and negative in #calc.round(r.neg_span.at(0))-#calc.round(r.neg_span.at(1)) s. Its absolute sign depends on the convention of the sign correction, but the change does not. Possible causes to check: trajectory (apogee), attitude, or drift direction of the structures.

- *The sign correction is supported by the data.* Without it, inside the significant blocks, segments after a positive spike give #n(r.pol_plus) and segments after a negative spike give #n(r.pol_minus): opposite signs, as expected from a sensor that flips.

- *Coupling follows wave activity (exploratory).* Per block, its size rises with E-field power at 10-40 Hz (Spearman $rho = #n(r.rho_e)$) and with density signal-to-noise ($rho = #n(r.rho_n)$). Before about 490 s the density power at 10-40 Hz equals its 60-75 Hz floor, so we cannot tell whether the coupling is absent there or only undetectable.

- *A caution on power-weighted coherence.* The classic coherence over the flight is #n(r.coh, d: 4), with $p = #n(r.p_shuffle, d: 3)$ against shuffled pairings but $p = #n(r.p_phase, d: 2)$ against a null that keeps the pairing and randomizes one phase per segment. Both signals strengthen together in the aurora, which inflates the estimate (effectively #calc.round(r.n_eff) of #r.segments segments carry it). We therefore rely on the equal-weight estimator.

- *Spectral slopes are a fragile test.* In 500-815 s, E-field slope minus density slope is #n(r.auto_diff_mean), 95% CI [#n(r.auto_diff_ci.at(0)), #n(r.auto_diff_ci.at(1))]; on the 2022 windows it ranges from #n(r.band_diff_min) to #n(r.band_diff_max) with the fit band. A working expectation of +2 assumes constant $T_e$, frozen-in structures and a measured component along the wave vector.

= 4. Next step: an automatic plasma-frequency trace

The density noise floor sets the limit of every result above, so the next step is to rebuild the trace from the raw HF spectra (about 50 million points for this flight, not in the public repository). In the existing trace, consecutive samples differ by 2 to 4 kHz (median) although the frequency bin is 0.3 kHz, so there is room to improve. Proposed method: (0) a robust outlier filter on the 2022 trace as the baseline to beat; (1) score each candidate frequency by the drop in power across it in the smoothed log spectrum; (2) find the best path through time with the Viterbi algorithm, below the upper-hybrid cutoff, with a jump penalty that is flat for small moves and grows only for large ones, so that it removes isolated spikes and dips without smoothing the fluctuations we want to measure; (3) return "no detection" where no cutoff stands out; (4) refine each point with a local step fit and attach an uncertainty. Validation: recovery of fluctuations injected into synthetic spectrograms, and a lower density noise floor at 10-40 Hz on real data; agreement with the hand-corrected trace (about 126,000 points) only catches gross errors, since that trace has its own biases. The coupling is computed only afterwards and is never used for tuning. This plan was written without access to the spectra and will need revision once they are inspected. The same pipeline would then apply to TRICE-2 Low, CAPER-2 and other rocket data sets.

= 5. Limitations

The instrument geometry (spin, probe configuration, calibration, frequency response) is inferred from the signal and needs confirmation; with a single channel the field cannot be despun. Blocks have a fixed number of segments, so their duration varies (28 to 109 s), and the block analysis was not planned in advance. The spike detector is conservative: it discards segments but never lets a spike in. Two original MATLAB scripts carry a "TRICE2-LOW" header; the data is believed to be from TRICE-2 High, to be confirmed.

#v(0.3em)
#text(size: 8.5pt)[*Reproducibility.* `make all` sets up Python 3.12, verifies data checksums, runs 14 tests (the exact MATLAB reproduction and synthetic checks of each method) and regenerates every number and figure in this note.]
