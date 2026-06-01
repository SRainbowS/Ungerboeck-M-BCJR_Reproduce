# Reduced-Complexity Receivers for Strongly Narrowband Intersymbol Interference Introduced by Faster-than-Nyquist Signaling

Adnan Prlja and John B. Anderson

Abstract—We propose new M-algorithm BCJR (M-BCJR) algorithms for low-complexity turbo equalization and apply them to severe intersymbol interference (ISI) introduced by faster than Nyquist signaling. These reduced-search detectors are evaluated in simple detection over the ISI channel and in iterative decoding of coded FTN transmissions. In the second case, accurate log likelihood ratios are essential and we introduce a 3-recursion M-BCJR that provides this. Focusing signal energy by a minimum phase conversion before the M-BCJR is also essential; we propose an improvement to this older idea. The new M-BCJRs are compared to reduced-trellis VA and BCJR benchmarks. The FTN signals carry 4–8 bits/Hz-s in a fixed spectrum, with severe ISI models as long as 32 taps. The combination of coded FTN and the reduced-complexity BCJR is an attractive narrowband coding method.

Index Terms—Faster-than-Nyquist, turbo equalization, minimum phase, intersymbol interference, BCJR, reduced complexity.

# I. INTRODUCTION

W E investigate the design and complexity of receiverswhen a convolutionally coded transmission is strongly band limited and the receiver is of the soft-input soft-output type. The modulation method is faster-than-Nyquist (FTN) signaling, i.e., linear modulation with a baseband pulse $h ( t )$ according to

$$
s (t) = \sqrt {E _ {s} / T} \sum_ {n} a _ {n} h (t - n \tau T), \quad \tau \leq 1 \tag {1}
$$

where $\left\{ a _ { n } \right\}$ are -ary independent and identically distributed Msymbols with zero mean and unit variance, $E _ { s }$ is the average modulation symbol energy, and $h ( t )$ is an arbitrary unit energy $T$ ( )-orthogonal pulse. The symbol time $\tau T$ , $\tau < 1$ , is shorter than $T$ 1; that is, the pulses are sent faster than allowed by the Nyquist orthogonality criterion and there is intersymbol interference (ISI) at the receiver. An additive white Gaussian noise (AWGN) channel with noise power spectral density $N _ { 0 } / 2$ follows $s ( t )$ . Form (1), with $\tau = 1$ , underlies most 2 ( )practical modulations.

The objective of this paper is two-fold: To explore iterative receivers for coded narrowband FTN signaling and to find new reduced-complexity BCJR algorithms for use in iterative

Paper approved by H. Leib, the Editor for Communication and Information Theory of the IEEE Communications Society. Manuscript received May 9, 2011; revised November 16, 2011 and March 5, 2012.

The authors are with the Dept. of Electrical and Information Tech. and the Strategic Center for High Speed Wireless Communication, Lund University, Lund, Sweden (e-mail: {adnan.prlja, anderson}@eit.lth.se).

Digital Object Identifier 10.1109/TCOMM.2012.070912.110296

![](images/5ffc7ed008c63c983928a0f19c3e917d1ce14458208bc256fe4b31d573cd5713.jpg)  
Fig. 1. Turbo equalization with a simple detection inner coder (dashed box).

decoding. The new algorithms follow the well-known Malgorithm idea, meaning that the BCJR recursions are based only on the $M$ dominant terms at each trellis stage. In an iterative decoder log likelihood ratio (LLR) values are passed around between two BCJR algorithms, and the quality of these strongly affects performance. Unfortunately, the M-algorithm degrades the LLR quality. However, signals can be preprocessed to make better use of the retained $M$ terms in the BCJR recursions. The contributions of this paper are improved minimum phase modeling and new M-BCJR algorithms that produce high quality LLRs, together with test results for coded FTN. The outcome is receivers for narrowband coding that work reasonably close to capacity with practical complexity.

Since FTN signals are continuous, a receiver contains a matched filter, sampler and possibly a post filter, which together reduce the signaling to a discrete-time convolution of the data $a _ { 0 } , a _ { 1 } , \dotsc$ (binary in this paper) with the ISI tap set $v = v _ { 0 } , v _ { 1 } , \ldots , v _ { \mathrm { m _ { T } } }$ . This provides a $2 ^ { m _ { T } }$ -state trellis for = 2the channel. Suitable receiver models are derived in Section II, which have the property that zero-mean IID Gaussians with variance $N _ { 0 } / 2$ are added to the discrete convolution 2values. We apply FTN signaling in two ways, by itself as an uncoded narrowband communication system and as the inner ISI mechanism in a coded system with iterative decoding. These are shown schematically in Fig. 1. The first embraces only the elements in the dashed box, and we will refer to it as “simple detection” of ISI. The second is often called turbo equalization [4].

Since the intentional ISI introduced by FTN with small

$\tau$ is severe, it is necessary to reduce the complexity of the ISI-BCJR block in Fig. 1. Coded FTN offers an attractive combination of bandwidth reduction and coding gain, but a reasonable receiver is needed. Reduced complexity can be achieved in two basic ways: By reducing the state size of the model $\pmb { v }$ , a reduced-trellis approach, or by reducing the search of a given trellis, a reduced-search approach. We will concentrate on the second and use the first as a performance benchmark. Early work with reduced search decoders primarily treats non-iterative applications where log likelihood ratios (LLRs) are not needed. A selection of papers on M- or other reduced complexity BCJRs is [6]–[11]. A factor graph based approach has been presented in [12]. In the remainder of this section we review some of the underlying ideas of FTN and reduced-search decoding.

In 1975 Mazo [1] pointed out that binary $\mathrm { s i n c } ( t / T )$ pulses sinc( )in (1) could be sent “faster” without loss of signal minimum distance. The demodulation asymptotic error probability is thus unaffected. He called this faster-than-Nyquist signaling. The pulse $h ( t )$ in our work is much narrower band than $1 / 2 \tau T$ ( ) 1 2Hz, and consequently severe ISI is introduced. Further, the interesting FTN cases have spectral zeros. FTN signaling has been extended in many ways: The modulation can be coded, it can be nonbinary, the pulse does not need to be · or sinc( )even orthogonal. The concept can be applied in frequency as well as in time [3], by placing OFDM-like subcarriers closer than orthogonality allows. An early study of receivers is [2]. In every case there will be a closest packing (a smallest $\tau$ and/or a closest subcarrier spacing) at which the signal minimum distance first falls below the isolated pulse value. This is called the Mazo limit to signaling with this $h ( t )$ and alphabet.

( )It has been known for some time that the receiver front end processing should provide the reduced-search detector with a minimum phase input (see [5], [7] and the more recent [11], [14]). A straightforward solution is to follow the matched filter/sampler by an all-pass filter that produces a max phase output, and then reverse the output frame. Min phase moves energy to the front of the ISI model, which directs the reduced search more efficiently (min phase will not improve a full $2 ^ { m _ { T } }$ -state Viterbi algorithm (VA) or BCJR). Energy focusing 2also aids reduced-trellis decoders, and in order to have a fair benchmark we will employ it there too. Section II discusses the min phase notion further, and proposes a novel extension to it that focuses energy more than simply calculating the mathematical min phase model.

An important role in this paper is played by the normalized minimum Euclidean distance $d _ { \mathrm { m i n } }$ between two signals of form (1). The error probability of maximum likelihood simple detection of the $\left\{ a _ { n } \right\}$ tends in logarithm to $\kappa Q ( \sqrt { d _ { \mathrm { m i n } } E _ { s } / N _ { 0 } } )$ , asymptotically in the ratio $E _ { s } / N _ { 0 }$ , where $\kappa$ ( )is a factor that depends on the most likely error events and whether bit error rate (BER) or event error rate (EER) is of interest. A union bound estimate formed from events at several distances near $d _ { \mathrm { m i n } }$ is a good estimate for simple detection of ISI at moderate to high $E _ { s } / N _ { 0 }$ . The details of bounding and distance finding appear in coding texts, e.g. ref. [22]. As the state size of a reduced-complexity VA or BCJR algorithm drops, its error rate will at some point depart from this ML estimate. Distance-based estimates are

thus essential for judging the minimum required size of an algorithm. For binary $a _ { n }$ , $d _ { \operatorname* { m i n } { } } \ \leq { } \ 2$ , and the special case $d _ { \operatorname* { m i n } } = 2$ 2means that the ISI can be detected with the same log = 2error probability as antipodal signaling except for the factor $\kappa$ ; that is, asymptotically the effect of the ISI can be removed. It also means that coded ISI can be decoded in a single turbo iteration; we will see a case of this in Sections V and VI.

In iterative decoding on the other hand, $d _ { \mathrm { m i n } }$ sets the error performance of the first ISI-BCJR pass. Above a certain $E _ { s } / N _ { 0 }$ , called the threshold, it can be shown that the iterations converge to the BER of the convolutional code alone over an antipodal (ISI free) channel with the same $E _ { s } / N _ { 0 }$ . Below this threshold, convergence is to a much higher BER.

The BCJR algorithm consists of forward and backward linear recursions, instead of the VA’s unidirectional add-compareselect. As such, its behavior is rather different from the VA’s. There is also a major difference between an algorithm that calculates full LLRs and one that makes decisions about bits, i.e., calculates the LLR sign. We find that accurate LLRs are essential in iterative decoding of FTN signals and that producing them is a considerable challenge for the Malgorithm. Earlier work on this subject appears in [11]. We find that the key to good quality LLRs is to add a third, low complexity recursion.

The algorithms in this paper apply to general ISI, but the intentional ISI introduced by FTN signaling is interesting for several reasons. First, it is severe, meaning that it has a combination of large state space, small $d _ { \mathrm { m i n } }$ , and $z$ -plane zeros on or near the unit circle. It is difficult to assess the effect of complexity reduction unless there is significant complexity to reduce. In Sections V and VI we will find that many of the ISI models in the literature provide a too easy target for reduction.

FTN is also interesting for theoretical reasons. The signals have a fixed power spectral density (PSD) shape, given by the Fourier transform of $h ( t )$ in (1). Such signals have a Shannon constrained capacity for the PSD, given by $\begin{array} { r } { \int ( 1 / 2 ) \log _ { 2 } [ 1 + } \end{array}$ $2 P | H ( f ) | ^ { 2 } / N _ { 0 } ] d f$ , where $P$ (1 2is the total power and $P | H ( f ) | ^ { 2 }$ 2 ( ) ] ( )is the PSD [17]. In general, this capacity cannot be reached by codes based on orthogonal pulses with PSD $P | H ( f ) | ^ { 2 }$ such as the coded modulations and turbo codes in common use. Studies of best convolutional codes have demonstrated that M-BCJR iterative decoders reach to 1–2 dB from this PSD capacity, with complexity and block length comparable to other iterative decoding [16].

A third reason for FTN signals is that they provide a proper experimental design for narrowband signaling. We explore the behavior of reduced BCJRs as the bandwidth drops and decoding becomes more complex. With narrowband signals receiver error performance is sensitive to the entire shape of the signal PSD, not just to a measure like 3 dB bandwidth. Minimum distance studies show that removal of only a small power from the outer spectrum can change the minimum distance of a signal set significantly; this is the “escaped distance” problem ([22], Chapter 6). These effects grow more pronounced as the bit density carried in bits/Hz-s grows. If a small extra power appears in the stopband—for example, through too-early truncation of the model $\textbf {  { v } }$ —receiver error rate can improve, and give a false test result for that model. FTN signals provide a way to increase the transmission bit

density, by reducing $\tau$ , while maintaining the same PSD shape.

Although any FTN pulse shape $h ( t )$ could be taken, in this paper $h ( t )$ ( )is the unit-energy root raised-cosine (rRC) $T$ - ( )orthogonal pulse with $30 \%$ excess bandwidth. Its spectrum is in fact zero outside $\pm 1 . 3 / 2 T$ Hz. Setting $\tau \ = \ 1$ gives 1 3 2the widely used rRC orthogonal pulse. As $\tau$ = 1drops below 1, pulses are sent “faster” but the PSD shape remains the same, namely, a raised cosine. The bit density is $2 / \tau$ data bits/Hz-s (taking 3 dB bandwidth as a scale unit). The asymptotic error rate remains $Q ( \sqrt { 2 E _ { s } / N _ { 0 } } )$ for $\tau \geq . 7 0 3$ , the Mazo limit. (Thereafter, it is $\approx Q ( \sqrt { d _ { \operatorname* { m i n } } ^ { 2 } E _ { s } / N _ { 0 } } )$ 703, where $d _ { \mathrm { m i n } }$ declines with $\tau$ ( ). The Mazo limit itself depends on the pulse excess bandwidth, i.e., $d _ { \mathrm { m i n } } ^ { 2 }$ will fall below 2 at a different $\tau$ Optimizing $h ( t )$ within a suitable framework is an interesting ( )future topic; some work of our own is [14].

The paper is organized as follows. Section II presents a suitable receiver front end and an improved discrete-time model, that yield white noise and easy control of spectrum and min phase. Section III reviews the BCJR. Section IV presents a benchmark offset VA and BCJR for severe ISI. Sections V and VI present and evaluate novel M-BCJR algorithms for simple detection and iterative decoding.

# II. DISCRETE-TIME SYSTEM MODELS

The conversion of continuous FTN signals to discrete time is considered in this section. Many methods are possible, and by choosing one we create the discrete-time signal model $\pmb { v }$ seen by the detector/decoder. In choosing a method for this paper, we have three priorities: Signals with spectral zero regions must be handled in an accurate, straightforward way, noise at the detector input should be white, and the model should be minimum phase.

We adopt the following model conversion to discrete time (henceforth called “conversion/model”), which was introduced in [15]. It assumes linear modulation by $h ( t )$ at rate $1 / \tau T$ and ( ) 1an AWGN channel, and then processes the signal according to:

Matched Filter Sample at nτT → Allpass filter $B ( z )$ Frame reverse Detector/Decoder.

The filter $B ( z )$ creates a max phase output, which is ( )reversed blockwise to form a min phase output. The matched filter is matched to some pulse $\phi ( t )$ and sampled at the faster rate $1 / \tau T$ . Let $\{ \phi ( t - k \tau T ) \}$ ), $k$ an integer, be an 1orthonormal basis for $h ( t )$ (, so that $\begin{array} { r } { h ( t ) = \sum c _ { k } \phi ( t - k \tau T ) } \end{array}$ where $\begin{array} { r } { c _ { k } = \int h ( t ) \phi ( t - k \tau T ) d t } \end{array}$ ( ) =. We choose $\phi ( t )$ (so that $\left\{ c _ { k } \right\}$ = ( ) ( )are the energy-normalized samples $h ( k \tau T )$ of $h ( t )$ . 1 Our $h$ is ( ) ( )infinite-response and time-symmetric, and there is a $K$ such that $c = \{ c _ { k } \}$ , $k = - K , \ldots , K$ will capture all but $\delta$ of the =pulse energy, any $\delta > 0$ . Since the $\{ \phi ( t - k \tau T ) \}$ are $\tau T .$ - 0 ( )orthonormal, the matched filter samples satisfy two important properties: Filtered noise samples are white Gaussian, and Euclidean distance between two noise-free continuous signals from (1) can be calculated from their samples.

Two other well-established conversion/models in the literature are the whitened matched filter (WMF) and the

Ungerboeck receivers. The WMF receiver filter is matched to $h ( t )$ and sampled each $\tau T$ . There follows a whitening ( )filter; its time-reversed output is min phase but it has stability problems when $h ( t )$ has spectral zeros. The WMF model is ( )applied in [14]. In the Ungerboeck receiver, the receive filter is again matched to $h ( t )$ and sampled each $\tau T$ , but there is ( )no whitening; a special processor works instead with colored noise. A BCJR of this type was explored in [24].

Turning to the allpass $B ( z )$ , we seek in the first instance the $B ( z )$ that makes $\pmb { v }$ ( )maximum phase. Allpass filters affect ( )neither the statistics of the noise (it is still white) nor the minimum distance of a signal set ([22], Chapter 6). This is true for any allpass. Max phase is achieved by a particular $B ( z )$ , the one that reflects outside the unit circle the zeros $\left\{ z _ { i } \right\}$ )of $\begin{array} { r } { C ( z ) = \sum c _ { k } z ^ { - k } } \end{array}$ that lie inside the circle; that is, the poles of $B ( z )$ =lie at $\left\{ z _ { i } \right\}$ and the zeros lie at $\{ 1 / z _ { i } \}$ . Zeros of $C ( z )$ ( )on the unit circle are not reflected.

( )With a reduced-search detector, there in fact exist $B ( z )$ that ( )improve the error rate even more than the mathematically correct $B ( z )$ . Reduced-complexity algorithms need a steep energy ( )growth in the model taps $\textbf {  { v } }$ . Suppose that $B ( z )$ produces a more rapid growth, but also a length- $K _ { p }$ low-energy precursor. Since the precursor energy is low, the algorithm can ignore it with almost no effect, i.e., it can work with a $\pmb { v }$ whose first $K _ { p }$ taps are set to zero. Consequently the detector is slightly mismatched to the true channel model. The key issue is: For a given reduced search does the better performance exceed the loss from the mismatch. This is a challenging optimization problem, and we lack space to pursue it here, but we list specific $B ( z )$ that we have found. An M-BCJR algorithm working with these $B ( z )$ can achieve the same bit error rate (with 2–4 times smaller $M$ . Such an improved model will be called super minimum phase.

The super min phase $B ( z )$ leads to significant improve-( )ment in the M-BCJR, but also to the reduced-trellis VA and BCJR benchmarks in the following sections. The physical decoder/detector remains the same, except that it runs $K _ { p }$ stages behind the present trellis stage and it computes branch labels ignoring the precursor. This $K _ { p }$ -delayed decoding is an essential element of the super min phase method.

We now give an illustration of min and super min phase models $\textbf {  { v } }$ for the $30 \%$ rRC FTN pulse stretched in time by $\tau = . 5$ . The central pulse samples are

$$
\begin{array}{l} \left\{c _ {k} \right\} = \left\{h (k \tau T) \right\} = \{\dots ,. 0 4 0, -. 1 0 9, -. 0 5 3, . 4 3 5, . 7 6 5, \\ . 4 3 5, -. 0 5 3, -. 1 0 9, . 0 4 0, \dots \} \tag {2} \\ \end{array}
$$

The max phase conversion of 61 of these, reversed, is plotted in Fig. 2. Now consider only the center-most 9 samples, the ones between the dots in (2). With only these, the reversed max phase conversion is

$$
\{. 3 7 5, . 7 4 2, . 5 0 0, -. 0 7 0, -. 2 1 6, . 0 1 4, . 0 7 7, -. 0 3 2, . 0 0 4 \} \tag {3}
$$

and the $B ( z )$ that creates it is

$$
B (z) = \frac {. 1 0 7 - . 5 6 1 z ^ {- 1} + z ^ {- 2}}{1 - . 5 6 1 z ^ {- 1} + . 1 0 7 z ^ {- 2}}. \tag {4}
$$

Even though the new model (3) energy rises faster it lacks the required PSD. If we instead filter the full model (2) by this $B ( z )$ , the outcome will have the correct spectrum, since

![](images/16d35d6004b07877fe0df37d30f7a7e45a15602f76e154877e281bd23fd44cdd.jpg)  
Fig. 2. Illustration of super minimum phase modeling at FTN $\tau = 1 / 2$ . Mathematical min phase response $_ { v }$ based on 61 pulse samples (circles); super min phase response (squares).

$B ( z )$ is an allpass. The significant parts of the outcome are ( )plotted (squares) in Fig. 2. The values at times $0 , \ldots , 8$ are 0 8nearly identical to (3); what is added is a precursor and the values at $9 , 1 0 , \ldots$ The latter points will not affect the M-9 10BCJR complexity. The precursor increases its complexity but if we can ignore it without damaging the error performance, this new $B ( z )$ will be a superior allpass. It will lead to a faster ( )rise of the main ISI model energy.

Figure 3 plots the improved FTN models $\textbf {  { v } }$ presented to the receiver processor for the main tests in this paper. The unit-energy models for $\tau = . 7 0 3 , . 5 , . 3 5 , . 2 5$ are respectively

$$
\begin{array}{l} \begin{array}{r c l} \boldsymbol {v} & = & [. 5 5 3,. 7 9 3, -. 0 8 4, -. 1 7 1, . 1 5 4, -. 0 6 4, . 0 0 6, . 0 1 0, -. 0 1 2, . 0 1 5, \\ & & -. 0 1 6, . 0 1 3, -. 0 0 8 ] \end{array} (5) \\ \begin{array}{r c l} \boldsymbol {v} & = & [ -. 0 0 5, -. 0 0 3, . 0 0 7, -. 0 1 1, -. 0 0 1, . 0 3 4, -. 0 1 9, . 0 0 3, . 3 7 5, . 7 4 1, \\ & & . 4 9 9, -. 0 7 0, -. 2 1 4, . 0 1 9, . 0 8 7, -. 0 2 0, -. 0 2 8, . 0 1 7 ] \end{array} (6) \\ \begin{array}{r c l} \boldsymbol {v} & = & [. 0 2 5, . 0 1 2, -. 0 2 4, . 0 0 8, . 1 9 1, . 4 6 4, . 6 2 3, . 5 0 6, -. 1 7 6, -. 1 2 3, \\ & & -. 1 9 6, -. 0 7 5, . 0 6 0, . 0 8 0, . 0 1 3, -. 0 3 5, -. 0 2 2 ] \end{array} (7) \\ \begin{array}{r c l} \boldsymbol {v} & = & \left[ -. 0 1 0, -. 0 1 3, -. 0 0 7, . 0 0 5, . 0 1 1, . 0 0 4, -. 0 0 8, . 0 0 1, . 0 6 0, \right. \\ & & \left. . 1 8 1, . 3 3 9, . 4 7 3, . 5 2 0, . 4 4 3, . 2 6 2, . 0 4 7, -. 1 2 0, -. 1 8 2, -. 1 3 8, \right. \\ & & \left. -. 0 3 7, . 0 5 5, . 0 9 2, . 0 7 0, . 0 1 8, -. 0 2 5, -. 0 3 7, -. 0 2 1, . 0 0 3, \right. \\ & & \left. . 0 1 6, . 0 1 2, . 0 0 0 4, -. 0 0 8 \right] \end{array} (8) \\ \end{array}
$$

The precursors are shown italic in (6)–(8); all detectors replace these with zeros and work at a delay $K _ { p }$ . The first $\tau$ is the Mazo limit for the $30 \%$ rRC $h ( t )$ . The last three ( )models are super min phase, with the allpass filter found from a search among $B ( z )$ obtained from truncations of $h ( t )$ . Note ( ) ( )that they have taps in the pattern [low energy precursor] $^ +$ [high energy part] $^ +$ [long decaying tail]. Insignificant taps before and after have been dropped.2 The test of whether too many taps have been dropped is the model spectrum, and these are plotted for each $\tau$ in Fig. 4. Compared to the ideal rRC spectra, spectral sidelobes must appear, but these are down at least 30 dB. Models with sidelobes down only 15–20 dB can have significantly better minimum distance than the true FTN

$^ 2 \mathrm { I n }$ tests the transmitted signal generation uses a few extra small taps, as insurance that the PSD is maintained.

![](images/3f50fb78d61687849fd8a9166279bb308e8e55ea9051a5e0d2f36e348a463aa1.jpg)

![](images/4c7d6a6e4186ce50b73b33106d50810537b257cdd1f84a1f1b60292018fe72ca.jpg)  
Fig. 3. Improved unit-energy discrete-time channel models, as seen by the VA or BCJR processor. FTN $\dot { \tau } = { . 7 0 3 , 1 } / { 2 } , { . 3 5 , 1 } / { 4 }$ .   
Fig. 4. Spectra of the channel models (dashed), compared to ideal $30 \%$ root RC spectra (solid); FTN $\tau = 1 / 4 , . 3 5 , 1 / 2 , . 7 0 3 , 1 .$ X-axis is $2 f T$ .

signals, and the receiver will show an artificially low bit error rate.

A detector branch label $\ell$ at trellis stage $n$ is generated from a $\pm 1$ data sequence $^ { a }$ by

$$
\ell = \sum_ {k = 0} ^ {m _ {T}} a _ {n - k} v _ {k} \tag {9}
$$

where the total memory $m _ { T }$ is the sum of the high energy and long tail lengths, and symbols $\ldots , a _ { n - 1 } , a _ { n }$ are the symbols following the precursor. The $\tau ~ = ~ . 5$ case generates mild ISI and a $50 \%$ = bandwidth reduction; $\tau = . 3 5$ is severe ISI and a reduction to $\approx 1 / 3$ = 35; the . case is very severe and a 25reduction to 1/4. The signal sets created by these ISIs have square minimum distances of 1.02, .56 and .20, which are energy losses of 2.9, 5.5 and 10.0 dB compared to antipodal signaling.

The following non-FTN discrete-time models from the literature will be used to compare with earlier work. They are much simpler, and the proposed M-BCJR will need to pursue only 2–3 paths to achieve near-ML performance. The model

$$
\boldsymbol {v} = \sqrt {. 4 5}, \sqrt {. 2 5}, \sqrt {. 1 5}, \sqrt {. 1}, \sqrt {. 0 5} \tag {10}
$$

features in early turbo equalization papers [4], [8]. It is min phase and $d _ { \operatorname* { m i n } } ^ { 2 } = 1 . 1 2$ ; the asymptotic VA equalizer EER is $\approx . 5 Q ( \sqrt { 1 . 1 2 E _ { s } / N _ { 0 } } )$ 2. The model

$$
\boldsymbol {v} = . 1 7 6 2, . 3 1 6 3, . 4 7 6 5, . 5 3 2 6, . 4 7 6 5, . 3 1 6 3, . 1 7 6 2 \tag {11}
$$

appears in several papers; it is min phas. . This tap set is said to have the least $d _ { \mathrm { m i n } } ^ { 2 }$ d has of an $d _ { \mathrm { m i n } } ^ { 2 } =$ $m = 6$ 2616set with binary input [13]. The model $\{ 1 , 0 , 1 , 2 , 1 , 0 , 1 \} / \sqrt { 8 }$ was studied in [11]. It has $d _ { \operatorname* { m i n } } ^ { 2 } = 2$ 1 0 1 2 1 0 1 8. In our receiver model it = 2would take its min phase form, which is

$$
\boldsymbol {v} = . 6 7 0,.. 3 6 6,.. 1 7 8,.. 4 4 3,.. 3 7 9. - . 1 0 2,.. 1 8 7 \tag {12}
$$

The super min phase idea is not useful with these short models.

# III. THE BCJR ALGORITHM

The BCJR algorithm computes probabilities of states and paths in a signal trellis, given the channel outputs ${ \textbf { 3 } } =$ $y _ { 1 } , \ldots , y _ { N }$ =and the apriori data probabilities. It efficiently calculates soft information in the form of LLRs, using the factorization

$$
P \left(s _ {n} = i, s _ {n + 1} = j, \boldsymbol {y}\right) = \alpha_ {n} [ i ] \Gamma_ {n} (i, j) \beta_ {n + 1} [ j ] \tag {13}
$$

where $\alpha _ { n } [ i ]$ is the forward trellis working variable for state $i$ at [ ]trellis depth $\mathscr { n } , \beta _ { \mathscr { n } + 1 } [ j ]$ is the backward variable for state $j$ at depth $n { \mathrel { + { 1 } } }$ , and $\Gamma _ { n } ( i , j )$ is the metric of the branch connecting states $i$ +1and $j$ Γ ( ). Starting from the initial all-zero state at the root of the trellis, the forward variable is computed recursively in a forward trellis pass according to

$$
\alpha_ {n + 1} [ j ] = \sum_ {i \in \mathcal {S}} \alpha_ {n} [ i ] \Gamma_ {n} (i, j) \tag {14}
$$

with the initialization $\pmb { \alpha } _ { 0 } = ( 1 , 0 , \dots , 0 )$ , where $s$ is the set = (1of states that can reach state $j$ 0)at stage $n + 1$ (in binary + 1transmission there are 2). Similarly, the backward recursion initialized with $\beta _ { N } = ( 1 , 0 , \ldots , 0 ) ^ { \prime }$ starts at the end of the = (1 0 0)trellis and proceeds to the root, computing at each depth

$$
\beta_ {n} [ i ] = \sum_ {j \in \mathcal {S}} \beta_ {n + 1} [ j ] \Gamma_ {n} (i, j). \tag {15}
$$

Now $s$ is the set of states reached from the stage- $n$ state $i$ The trellis branch metric is

$$
\Gamma_ {n} (i, j) = \left[ P \left(a ^ {\prime}\right) / \sqrt {\pi N _ {0}} \right] \exp \left[ - \left(y _ {n} - \ell_ {i, j}\right) ^ {2} / N _ {0} \right] \tag {16}
$$

where $\ell _ { i , j }$ is the label (9) on the branch from state $i$ to $j$ and $P ( { a } ^ { \prime } )$ is the apriori probability of the symbol $a ^ { \prime }$ that causes ( )the transition. In a reduced recursion, many of the smaller $\alpha _ { n }$ or $\beta _ { n }$ components are set to 0.

The product of the $\left\{ \alpha _ { n } \right\}$ and $\{ \beta _ { n } \}$ produce the set $\{ \lambda _ { n } \}$ through $\lambda _ { n } [ j ] = \alpha _ { n } [ j ] \beta _ { n } [ j ]$ , $j$ a state at stage $n$ . From these [ ] =we obtain LLRs via

$$
\operatorname {L L R} \left(a _ {n}\right) \triangleq \ln \frac {P \left[ a _ {n} = + 1 \mid \boldsymbol {y} , \mathrm {L L R} _ {\mathrm {i n}} \right]}{P \left[ a _ {n} = - 1 \mid \boldsymbol {y} , \mathrm {L L R} _ {\mathrm {i n}} \right]} = \ln \frac {\sum_ {j \in \mathcal {L} _ {+ 1}} \lambda_ {n} [ j ]}{\sum_ {j \in \mathcal {L} _ {- 1}} \lambda_ {n} [ j ]} \tag {17}
$$

Here $\mathcal { L } _ { \pm 1 }$ are the sets of states reached by $a _ { n } ~ = ~ \pm 1$ , for which nonzero $\alpha$ and $\beta$ = 1have both been found. A problem in a heavily reduced calculation is that one or both of $\mathcal { L } _ { \pm 1 }$ are often empty. Both being empty at a reasonable complexity never occurred in our tests. If only one set is empty, the numerator or denominator of (17) must be replaced by some backup method. This is a characteristic of reduced-search BCJRs and in Section V we propose a novel method for the M-BCJR which considerably improves coded FTN error performance.

# IV. REDUCED-TRELLIS BENCHMARKS: THE OFFSET BCJR AND VITERBI ALGORITHMS

This section sets up reduced-trellis benchmark detectors based on the full VA or BCJR, applied here to simple detection of uncoded ISI. In Section VI we will compare the benchmark error performance to the proposed M-BCJRs. Finding a fair benchmark is challenging. In fact, considerable trellis reduction is possible without significant error rate loss. Since algorithms that process reduced trellises are quite simple, the M-BCJR state search needs to be small in order to compete. We also want to distinguish the two types of complexity reduction, and see how they compare in FTN. Further requirements for a benchmark are that it must work within the constraints of Section II, which are white noise and error performance implied by the full signal set $d _ { \mathrm { m i n } }$ .

A key to reducing the state space of the VA or BCJR is to favor high-energy model taps, if it can be done simply, without increasing the error rate. We assume that the algorithm is preceded by the Section II conversion/model, so that $\textbf {  { v } }$ is min/super min phase, with energy focused near the present symbol. The following offset receiver will then reduce the complexity induced by the low-energy tail: Instead of generating trellis branch labels as in (9), form them instead from

$$
\ell = \sum_ {k = 0} ^ {m} a _ {n - k} v _ {k} + \sum_ {k = m + 1} ^ {m _ {T}} a _ {n - k} v _ {k}. \tag {18}
$$

Symbols $a _ { n - m } , \ldots , a _ { n - 1 }$ comprise the size- $m$ reduced VA/BCJR main state, and stem from high energy symbols, while the $a _ { n - m _ { T } } , . \cdot \cdot , a _ { n - m - 1 }$ comprise the offset state. The second term is an offset to the label $\ell$ created by early symbol history. A set of offset symbols can be associated with each main state but its symbols do not form part of the algorithm’s state variable. However, all $m _ { T } { + 1 }$ taps contribute to a label. In +1the add-compare-select step of the benchmark offset VA, the offset states of the survivors, together with the oldest main state bit, become the offset states for each new main state. Trellis searching focuses on high energy taps while small taps contribute only to the labels.

This sort of trellis reduction was devised in the 1970s [18] as a way to handle large state spaces of long-response systems, and was applied to ISI problems by several authors in the 1980s. In the best known [19], Duel-Hallen and Heegard calculate $d _ { \mathrm { m i n } }$ for the VA receiver as a function of the main state size $m$ . Studies of the VA then and now [14] show that under narrowband ISI a large truncation is possible without significant loss in $d _ { \mathrm { m i n } }$ . Offset BCJR receivers have been studied since the mid 1990s, although not for narrowband ISI. A major work is Colavolpe et al. [8], which gives a full list of references.

A different strategy to reduce trellis size is to add nonallpass prefiltering to the conversion/model. Even though they appear promising [8], [20], [21], these methods color the noise and reduce $d _ { \mathrm { m i n } }$ , and so we do not explore them further. Ref. [20] presents error rate results, to which we compare in the sequel.

Next we explore the benchmark VA performance. The offset VA is quite different from the BCJR and its benchmark

![](images/edebe42726198757811a27893a19add04d0ca7200fc4b024853d64a40b855e9c.jpg)  
Fig. 5. Error event rates versus $E _ { s } / N _ { 0 }$ for the offset VA with mathematical (dots) and super (solid) min phase channel models. FTN signals with $\tau =$ 1/2. Main state memory 2, 4, 6.

achieves near-optimal error performance with memory $m \ : 1 -$ 2 smaller under our severe ISI. We first consider the $\tau = . 5$ = 5FTN case in Figure 5. The offset VA is the standard kind [19], which associates an offset state with each main state.3

We wish to find the smallest $m$ that leads to essentially ML performance. The solid curves plot EER for the super min phase model (6) at main state memories 2,4,6; it is clear that $2 ^ { 4 } { - } 2 ^ { 6 }$ states are needed, and so this benchmark state size 2 2is about 32. The dotted curves show the offset VA with the mathematical min phase model derived from (2) at the same $m$ . Clearly this modeling is worse than the super min phase, especially for short memories.

Now we consider the BCJR benchmark algorithm. We apply recursions (14)–(15) to the main state in (18), computing the $2 ^ { m + 1 }$ branch labels with the help of the second-term label 2offsets (similarly to [8]). With long, narrowband ISIs of the FTN type we find that certain changes to the offset procedure improve performance, and in the interests of a fair benchmark comparison, we describe them briefly. We find that only a single offset should be associated with all the main states, not a different one for each state, as in the standard offset VA. Furthermore, the symbols used to compute the offset in (18) should be soft values, not $\pm 1$ . A solution to this comes from the definition of $\alpha$ 1in the BCJR, which has the form $\alpha _ { n } [ j ] \triangleq$ $P$ Observe $y _ { 1 } , \ldots , y _ { n } \cap \mathrm { I S I }$ is state $j$ at $n ]$ . Summing $\alpha _ { n }$ [ ]over [the states in $\mathcal { L } _ { + 1 } ^ { n - m }$ ] , the set of states that have entering symbol $+ 1$ at stage $n - m$ , we get

$$
\begin{array}{l} \pi_ {+ 1} = \sum_ {j \in \mathcal {L} _ {+ 1} ^ {n - m}} \alpha_ {n} [ j ] \\ = P \left[ \text {O b s e r v e} y _ {1}, \dots , y _ {n} \cap + 1 \text {s e n t a t} n - m \right] \tag {19} \\ \end{array}
$$

and similarly for $\pi _ { - 1 }$ . The probability of $+ 1$ at the oldest main state stage is then estimated as $\hat { p } _ { + 1 } = \pi _ { + 1 } / ( \pi _ { + 1 } + \pi _ { - 1 } )$

3The test setup is: Size 800 frames of random $\pm 1$ data, with enough frames to give 40–100 error events; frames are terminated before and after by $m _ { T }$ $_ { \cdot + 1 } \cdot$ symbols. The VA output symbol is decided $m _ { T } + 3 5$ symbols before the present trellis stage. Error events are taken to begin when the receiver output state splits from the transmitter state path and to end after 5 output data are correct. BER is about 3 times EER at higher $E _ { s } / N _ { 0 }$ and 4–5 times at lower.

![](images/8d153edca858798b280c010d1e2ad7e61288f5017070a1536b50327bc8d2f092.jpg)  
Fig. 6. Benchmark error event rates vs. $E _ { s } / N _ { 0 }$ for simple ISI detection with offset VA (solid) and single offset BCJR (dash dot) at main state memory $m$ . Heavy lines are Q-function estimates.

and similarly for $\hat { p } _ { - 1 }$ . These enable early decisions about $a _ { n - m }$ ˆ, and the single offset term can be built up from them. Although not as reliable as those based on the tworecursion BCJR, they are good enough for calculating a single label contribution from small taps. Furthermore a simple and effective soft decision about $a _ { n - m }$ is its expected value, which is $\hat { a } _ { n - m } = \hat { p } _ { + 1 } - \hat { p } _ { - 1 }$ . A highly likely $\pm 1$ is respectively $\approx \pm 1$ = ˆ ˆ 1and a completely uncertain symbol is 0, meaning that 1the tap is ignored in the offset computation. The single soft offset innovation can improve the BER of an offset BCJR used for simple detection by 10 fold. Some further details and alternative BCJRs are given in [15], but the best performing one is the single soft offset BCJR.

Final symbol decisions are made with both $\alpha$ and $\beta$ . From the $\alpha$ recursion, the benchmark single offset BCJR obtains a tentative estimate of $\textbf { \em a }$ , which can be used to compute offsets in the $\beta$ recursion. Note that in the $\beta$ recursion the channel response is backwards and the long tail symbols are not yet reached.

The general behavior of this single offset BCJR in simple detection, as a function of the offset and precursor size, differs little from the offset VA, although the BCJR requires 1–2 extra units of main state memory to achieve the same error performance. Figure 6 compares benchmark VA and BCJR EERs at main state $m = 2 , 4 , 6 , 7$ , for uncoded FTN with $\tau = . 5 , . 3 5$ = 2 4 6 7and super min phase models (6)–(7). The same test = 5 35setup as in Fig. 5 is adopted and both form tentative symbol estimates at a delay $m$ , using (19). The heavy lines are the Q-function estimates, based on the full model $\pmb { v }$ .4

Our study shows that the single offset BCJR needs only 32 states $( m \ : = \ : 5 )$ ) at $\tau ~ = ~ . 5$ and 64–128 $m = 6 { - } 7 ,$ ) at $\tau = . 3 5$ = 5 = 5 =. The offset VA needs somewhat less.5 With $\tau = . 2 5$ = 35(not shown) the offset VA needs about $2 ^ { 9 }$ = 25states, and the BCJR

$^ 4 \mathrm { A }$ distance study shows that the $d _ { \mathrm { m i n } }$ -causing error difference sequence is $2 , - 2 , 2$ for $\tau = . 5$ , with coefficient 1/4 (see [14] for details, and [22] for a general treatment). Thus the full-state VA has EER ≈ . $2 5 Q ( \sqrt { 1 . 0 2 E _ { s } / N _ { 0 } } )$ . The $\tau = . 3 5$ and .25 cases are more complex. The most probable error events at $\tau = . 3 5$ have differences $2 , - 2$ and $2 , - 2 , 0 , 2 , - 2$ ; these combine to yield an EER close to . $3 5 Q ( \sqrt { . 5 6 E _ { s } / N _ { 0 } } )$ for ${ E _ { s } } / { N _ { 0 } } = 1 0 \mathrm { { - } } 1 5$ dB.   
5This $m$ is also predicted for the VA by the reduced-trellis $d _ { \mathrm { m i n } }$ algorithm in [19].

about $2 ^ { 1 3 }$ . These numbers are a benchmark for the M-BCJR 2results in Section V. Observations for the non-FTN ISI tap sets (10)–(12) are given at the end of Section V.

# V. PROPOSED M-BCJR ALGORITHMS AND THEIR SIMPLE DETECTION PERFORMANCE

In this section we propose three new M-BCJR algorithms and test them in simple detection of ISI. The basic Malgorithm for reduced-search of trees and trellises is well known. As a general procedure, the algorithm proceeds breadth-first through a tree structure of metric values, keeping only the dominant $M$ paths at each tree stage. In a BCJR, it is applied once each to find the dominant $M \ \alpha _ { n } [ i ]$ and $\beta _ { n } [ j ]$ , near the values that a full BCJR would find at $n$ [ ] [ ]. For moderate to strong ISI, the branch metric matrices $\Gamma _ { n }$ are very sparse, and most non-zero elements are very small. A useful view is that the M-search implements a sparse matrix calculation in which the vector $\alpha _ { n }$ or $\beta _ { n }$ at each stage is limited to $M$ active components.

The most straightforward M-BCJR now follows. It works well in simple detection of ISI, and hence we call it the simple detection M-BCJR. Recursions start and end in state 0 (all $+ 1$ symbols). Inputs to the algorithm are the noisy channel outputs $\textbf {  { y } }$ and apriori probabilities of the symbols. Outputs are the signed LLR values in (17). The list of $M$ dominant paths consists of two sublists, one containing $\alpha$ or $\beta$ values at stage $n$ and one containing the corresponding trellis states. It is straightforward to extend the algorithm to non-binary alphabets.

Forward Recursion for $\alpha$ . Starting at $n = 0$ , perform at stage $1 , 2 , \ldots , N - 1$ :

1 21) The $\alpha$ 1recursion in (14) is computed from the $M$ nonzero values retained in $\alpha _ { n }$ . There are $M$ outcomes corresponding to symbol $a _ { n + 1 } = 1$ and $M$ to $- 1$ ; only the $2 M$ corresponding $\Gamma _ { n }$ = 1elements are computed.   
2) Trellis paths may merge at stage $n + 1$ . Merges are + 1detected and removed, leaving one survivor whose $\alpha$ value is the sum of the two incoming values.   
3) The best $M$ of the remaining paths are found. These are stored for the next stage and for the $\beta$ recursion.

Backward Recursion for $\beta$ . Starting at $n = N$ , the end of the channel block, perform at stage $N , N - 1 , \ldots , 2$ :

4) The $\beta$ 1recursion in (15) is computed from the $M$ nonzero values retained in $\beta _ { n + 1 }$ . There are $M$ outcomes corresponding to symbol $+ 1$ and $M$ to $- 1$ ; only the $2 M$ corresponding $\Gamma _ { n }$ elements are computed.   
5) Trellis paths may merge at stage $n$ . Merges are detected and removed, leaving one survivor whose $\beta$ value is the sum of the two incoming values.   
6) The best $M$ of the remaining paths are found, subject to the following condition: $\beta$ paths must be kept if their state and stage overlap with that of a stored $\alpha$ . The M-list is then completed with non-overlapping paths.   
7) Compute the LLR from (17). If $\mathcal { L } _ { + 1 }$ or $\mathcal { L } _ { - 1 }$ is empty, the respective $\lambda$ -sum in (17) is set to $\epsilon$ , where $\epsilon$ is a reserve value set in advance.

We can find no evidence that the offset state idea is needed in the M-BCJR; it should simply retain all $m _ { T }$ state symbols for each of the $M$ paths. It is essential, however, that the M-BCJR ignores the precursor symbols, and it is therefore slightly mismatched to the channel. The overlap in step 6 never failed to occur in our tests. The merges removed in steps 2 and 5 are unlikely, and these steps can be removed without significant performance loss. The idea of a reserve $\epsilon$ in step 7 and of pursuing $\beta$ paths that overlap $\alpha$ paths in step 6 were both proposed in [11] (the $\alpha$ path list was called the “survivor list”). However, we give the overlapping paths only first priority; other $\beta$ paths are extended for a total of $M$ . The efficiency of this strategy may be seen by observing the search dynamics. During most of the transmission, there are only 1–2 paths in the search overlap, and almost always one of these is correct. Errors occur during rare noise bursts, but it is precisely here that the search is chaotic; the extra $\beta$ paths are needed in case they merge to $\alpha$ paths a few stages later.

Comments on complexity. We comment only on the approximate order of the computations. Steps 3 and 6, which find the best $M$ , are equivalent to finding the median of a group of items. An important property of median finding is that its computation is linear in $M$ . The M-algorithm does not order a list, which would require order $M \log M$ . In keeping with this, logwe take a true M-algorithm to be one where all computation is of order $M$ . The search for the median is thus implemented in order $M$ , but so is also removal of state merges in steps 2 and 5 and finding the overlap of $\alpha$ and $\beta$ in step 6. The key to the last two is keeping all path lists in state order, which is itself a linear operation. We omit the details of these linear procedures. Since there are two recursions,6 computation has the approximate order $2 M$ , with $M$ the number of trellis states 2visited at each stage; by this measure the reduced BCJR has twice the complexity of the reduced VA.

Figure 7 plots the EER for this M-BCJR algorithm used as a simple detector at the ISI intensities τ . , . , . . = 5 35 25The algorithm decides symbols from the LLR sign at its output. Heavy lines show Q-function estimates. For comparison, performances are plotted for the 256- and 4096-state offset VA, for $\tau = . 3 5$ and . respectively. The simple detection = 35 25M-BCJR can perform better than the offset VA, especially at $\tau = . 2 5$ , because a practical VA cannot be large enough to deal = 25with every detail of the severe ISI. The M-BCJR needs only $M = 3 , 7 , 2 0$ respectively for the three FTN cases at higher $E _ { s } / N _ { 0 }$ 7 20, and somewhat more at lower. The appearance of such an upper limit to $M$ is typical of M-algorithm searching of code and modulation trellises. Not shown in the figure are results for $\tau = . 7 0 3$ , the Mazo limit; here only $M = 3$ is required.

Comparisons to Earlier Work. For mild ISI, some comparisons are possible with the literature. Magarini et al. [20] investigate alternatives to the benchmark offset VA, using ISI model (11). We obtain the same EER and BER results as they do for the offset VA (Figs. 6 and 7 of [20]), which here needs $m = 4$ (16 states); they are able to improve this benchmark

![](images/70b939c64e365600007d36334286917011bfcfbba767e913b211cd358f8e543b.jpg)  
Fig. 7. Error event rates for simple ISI detection vs. $E _ { s } / N _ { 0 }$ in dB for the simple detection M-BCJR (dotted lines); shown for comparison are offset VA (dash-dot) and Q-function estimate (solid).

performance somewhat with a non-allpass prefilter receiver. However, the simple detection M-BCJR with only $M = 7$ = 7improves upon their prefilter result by 1.5 dB and in fact lies on the full ML bound for BER given in [20].7

Fertonani et al. [11] consider ISI (12) in a turbo equalization configuration, but $d _ { \mathrm { m i n } } ^ { 2 }$ is 2 and as discussed in Section I iterative decoding is in principle not needed. The M-BCJR with $M = 6$ applied to simple detection of ISI (12) achieves = 6BER close to $\ 3 Q ( \sqrt { 2 E _ { s } / N _ { 0 } } )$ , which is the asymptotic esti-3 ( 2mate from distance analysis.

When an accurate LLR is needed—as in iterative decoding in the next section—the simple detection M-BCJR is not sufficient. A serious problem is that an empty $\mathcal { L } _ { \pm 1 }$ set normally occurs when $E _ { s } / N _ { 0 }$ takes practical values, since then the M-search is quite sure of the correct symbol; there is no estimate of the LLR magnitude at all, other than the  set apriori. Several solutions exist in the literature, and we and our colleagues have tested many more. The only practical one that we have found adds a third, low complexity recursion, whose purpose is to produce a backup LLR magnitude when the two M-recursions do not. We propose the following, called the backup M-BCJR. It replaces step 7 with:

New step 7. Decide the symbols from the sign of (17), noting when $\mathcal { L } _ { \pm 1 }$ is empty. In a third recursion, compute a symbol probability from the αs only, as follows. From each node of the decided symbol path, trace forward through the ISI trellis a certain length of stages; αs that stem from the decided branch form the probability of one symbol outcome and αs in the “incorrect subset” of the node form the probability of the other outcome. The traces are performed with a small Msearch of size $M _ { B }$ (typically $M _ { B } = 2$ works well).

= 2The necessary search for all the decided nodes at once can be arranged in a simple way. The search gives a backup estimate of $P [ a _ { n } = + 1 ] / P [ a _ { n } = - 1 ]$ , to be used when one (or both) of $\mathcal { L } _ { \pm 1 }$ [ = +1] [ = 1]is empty; otherwise (17) is used. A sketch of the entire backup M-BCJR procedure with $M = 3$ is illustrated in Fig. 8. First, the $\alpha$ = 3recursion is performed, then the $\beta$ . The $\beta$ paths, shown dotted, must follow the $\alpha$ paths as a first priority, and in this case they are shown as overlapping. Shown

![](images/85322ffc81ddd6921bb55d1f876299c16d15f839065f9a00a7cabc9860b5d610.jpg)  
Fig. 8. Backup M-BCJR example, showing $_ \alpha$ and $\beta$ recursions, hard decision path, and backup recursion.

third is the decided path that results from the whole first two recursions. Finally, a backup search with $M _ { B } = 2$ is shown for a branch that was decided to be $- 1$ .

1The backup LLR values can be noisy in a heavily reduced $M$ -search, especially in the early iterations of a turbo decoder, and therefore a useful technique is to smooth them. A simple moving average filter, such as . $2 z + . 6 + . 2 z ^ { - 1 }$ , can improve 2 + 6 + 2the BER of the iterative decoder, if only applied to the backup values in the first iteration. This third scheme is called the smoothed backup M-BCJR.

# VI. TURBO EQUALIZATION

Now we evaluate the BER performance of turbo equalization when the smoothed backup M-BCJR performs the ISI detection. Whereas only the sign of the LLRs was needed for simple ISI detection, turbo decoding requires reasonably accurate values, especially in the early iterations. This is provided by the backup M-BCJR. We will compare it with two reduced-trellis benchmark BCJRs, the memory- $\mathbf { \nabla } m$ single offset

![](images/b75fddc9bd0858cf3722f25241d9c9e58a7854a2754060743bedeba925db8888.jpg)  
Fig. 9. Turbo equalizer BER vs. $E _ { b } / N _ { 0 }$ for $\tau = 1 / 2$ , comparing single offset BCJR (dashed), smoothed backup M-BCJR (solid) and truncated BCJR (dotted) for different complexities.

BCJR from Section IV and a truncated BCJR that simply calculates its branch labels based on the $m _ { \mathrm { t r } } + 1$ dominant + 1taps with no offset to the labels. EXIT charts [23] are used to monitor convergence behavior. Low-quality LLRs affect the stability and convergence of the iterative detector and some way needs to be found to keep it under control when the complexity is heavily reduced. As a complement to the backup M-BCJR, it is beneficial to scale the likelihoods passed around the turbo loop by a “gain” $g \le 1$ (gains were suggested in 1[11]). We scale the extrinsic LLRs by $\sqrt { g }$ before each BCJR.

The turbo equalization setup is as follows: A block of $N$ information bits is encoded by the (7,5) rate 1/2 feed-forward convolutional code, generating a length $2 N$ codeword. The encoded sequence feeds a size $2 N$ 2random interleaver whose 2output is mapped to binary symbols $( 0 ~  ~ + 1 , 1 ~  ~ - 1 )$ . 0 +1 1 1The signal is terminated so that the transmission begins and ends in the all $+ 1$ ISI state. Iterative decoding as in +1Fig. 1 is performed, applying one of the three BCJRs as the ISI equalizer. All three ignore precursors when forming labels. In the smoothed backup M-BCJR smoothing is applied only at the first iteration with the smoother $\{ 1 , 3 , 1 \} / 5$ . The 1 3 1 5convolutional decoder is a full-state BCJR (4 states). In this section the signal-to-noise ratio (SNR) is defined as $E _ { b } / N _ { 0 }$ where $E _ { b } = 2 E _ { s }$ .

The component decoders exchange soft information in the form of LLRs, hopefully converging to a decision about the data. A complete loop is an “iteration”. The block length is $N = 1 2 0 0 0$ and 20 iterations are performed (60 for $\tau = . 2 5$ ). = 12000Fewer iterations and shorter blocks $\langle N \approx 1 0 0 0 \rangle$ = 25 are more 1000practical in hardware, and these performed almost as well, but more care is needed to assure loop stability. Decoder tests are run until $\geq 5 0$ blocks with errors occur.

The values $M$ in our tests are chosen with practical receivers and good performance at reasonable SNR in mind. Measurements show that one of the sets $\mathcal { L } _ { \pm 1 }$ is empty at $5 0 \mathrm { - } 8 0 \%$ o f the trellis stages and consequently some sort of LLR reserve procedure is essential. Higher quality LLRs and stabilizing loop gains allow smaller $M$ . Our test show that the best loop gains lie near . for $\tau = . 5$ and .35, and . for $\tau = . 2 5$ .

![](images/0f88c3fd36b4f4957cd065c0e79f59be2b13dbe3b988a412c2a2181fb73ea7b2.jpg)  
Fig. 10. Turbo equalizer BER vs. $E _ { b } / N _ { 0 }$ for $\tau = 0 . 3 5$ ; single offset BCJR, smoothed backup M-BCJR and truncated BCJR as in Fig. 9.

We use these values. They are much larger than those in [11]. The ultimate turbo goal, since no precoding is employed, is to reach the no-ISI performance of the convolutional code, which is shown in the plots as a heavy dashed ‘CC’ line.

A constant $M$ is employed in this paper. However, the first few iterations are the most important, and both $M$ and $g$ should vary with the iterations; this will be explored in a future paper. The first iteration is precisely the simple detection of Section V, so a suggestion for at least the starting $M \mathrm { s }$ are the values found there.

Figure 9 shows BER results for the smoothed backup M-BCJR, single offset and truncated BCJRs when $\tau = . 5$ with = 5taps (6). The relatively mild ISI is not difficult for the first two, but the truncated BCJR suffers from energy loss caused by a too-early truncation and fails to converge to the CC line at high SNR. For $M \ = \ 2$ the smoothed backup M-BCJR performs better than the memory-1 single offset BCJR, which has similar complexity, and it achieves very nearly the CC-line BER at $\mathrm { S N R } \approx 5$ dB. For higher complexities the smoothed 5backup M-BCJR and single offset BCJR are similar but the last is clearly superior to the truncated BCJR, which shows that some of the long-tail taps cannot be ignored.

The situation changes when the FTN signaling rate increases. Figure 10 plots the $\tau = . 3 5$ case, which is a much = 35more severe ISI. The smoothed backup M-BCJR efficiently removes the ISI even when $M ~ \leq ~ 5$ . It is superior to the 5single offset BCJR for all complexities. The reduced trellis of the truncated BCJR is now much smaller than the effective state space of the ISI, causing severely degraded BER. Even for $m _ { \mathrm { t r } } = 5$ the truncated BCJR is unable to eliminate the = 5effects of the intense ISI. Accounting for the long tail taps makes it possible for the single offset BCJR to achieve the CC performance even though the number of main states is equivalent to that of the truncated BCJR. The reduced-search M-BCJR clearly prevails at this higher ISI intensity. Its turbo convergence threshold can be determined through a study of EXIT charts, with the system converging to the CC line when there is an open tunnel between the EXIT curves. Figure 11 shows the case $\tau = . 3 5$ , $M = 5$ and $E _ { b } / N _ { 0 } = 7$ dB, for = 35 = 5 = 7which there is a narrow tunnel; Figure 9 verifies that there

![](images/c5ae88136aaddf27a9bf19a2eb86b382d0ec83dee1353fcb5f385046311831b8.jpg)  
Fig. 11. An EXIT chart at $E _ { b } / N _ { 0 } = 7 ~ \mathrm { d B }$ , showing extrinsic vs. apriori information for block length 12000. Dashed curve represents the smoothed backup M-BCJR with $M = 5$ and $\tau = 0 . 3 5$ ; solid curve shows the (7,5) outer convolutional code.

![](images/65be7d5b0a60202c891e42e6ae34e538fe3e1a7a61171de56ff8f3578121e103.jpg)  
Fig. 12. Turbo equalizer BER vs. $E _ { b } / N _ { 0 }$ for $\tau = 1 / 4$ ; smoothed backup M-BCJR with $M = 2 0 – 1 0 0$ . CC reference lies far below the SNR axis.

![](images/de5d843c1f8bf7db669cd05d76767578d699bf491d98546cd99994d9745d7eb8.jpg)  
Fig. 13. Turbo equalizer BER vs. $E _ { b } / N _ { 0 }$ for backup M-BCJR with $M _ { B } =$ $0 , 1 , 2 , 4$ , with comparison to algorithms from [9] and [11]. The second is tested with both math and super min phase model (7). FTN $\tau = . 3 5$ and $M = 6$ .

indeed is convergence to the CC line with $M = 5$ for the first time at about 7 dB.

Turbo equalization (60 iterations) for the extreme ISI case with $\tau = . 2 5$ and the 32-tap channel model (8) is shown in = 25Fig. 12 for several $M$ . Here the stability of the turbo loop is very sensitive to the gain $g$ and the block length. Lengths less than 12000, smaller $M$ and $g$ too large severely degrade the BER. The $M$ in the smoothed backup M-BCJR needs to be in the range 25–100, compared to 20 in simple detection of ISI and many thousands for the two reduced-trellis benchmark BCJRs. Depending on $M$ , sharp convergence thresholds lie in the range 9–10.5 dB; after the last point shown for each curve there is a sudden drop to the CC line, which is located far below the SNR axis.

As further insight into the role of the backup recursion, Fig. 13 shows the M-BCJR BER for several backup $M _ { B }$ when $\tau = . 3 5$ . A heavily reduced search with $M = 6$ leads to empty $\mathcal { L }$ = 35sets for $80 \%$ of the LLRs. The case $M _ { B } = 0$ corresponds = 0to the simple detection M-BCJR with a small, fixed reserve value $\epsilon$ whenever an $\mathcal { L }$ set is empty. Its BER performance is 3 dB worse than the backup M-BCJR with $M _ { B } = 4$ . However, = 4most of the performance gain is obtained with only $M _ { B } = 2$ . = 2We also compare to published M-BCJRs of which the most important one appears in [11] and is similar to our simple detection M-BCJR without steps 4 and 6. Its BER performance is shown with both mathematical min and super min phase models.8 The figure shows that both the backup recursion and the super min phase ideas are needed and that the gain from both innovations is about 4 dB. A final comparison in Fig. 13 is made to the $\mathbf { M } ^ { * }$ -BCJR algorithm of Sikora and Costello [9]. It has similar performance to the backup M-BCJR with $M _ { B } = 2$ , but it is much more complex.

= 2Comparisons to Earlier Work. Colavolpe et al. [8] report results for turbo equalization with tap set (10) and the 16-state recursive systematic convolutional code (23,35). They use an offset BCJR, block size 2048 data bits, and loop gain $g = . 1 5$ . = 15Their system needs 6 iterations and a 16-state BCJR to reach the CC line BER at 5 dB, although 8 states performs well (Fig. 5 of [8]).9 With the same setup (except $g = . 4 4 )$ ), we find that the smoothed backup M-BCJR with $M = 3$ 4needs = 3only 4 iterations at SNR 5 dB and 3 iterations at 6–7 dB.

Fertonani et al. [11] test several M-BCJR decoders with ISI (12) and our feed-forward (7,5) convolutional code, as mentioned in Section V. Their M-BCJRs need 6–8 paths and 20 iterations to reach the CC line (see Fig. 5, [11]). We find that the smoothed backup M-BCJR needs $M = 3$ and only $5 -$ = 311 iterations, depending on the SNR. The poorer performance in [11] probably stems from the lack of a backup recursion and any min phase conversion.

# VII. CONCLUSION

We have proposed and investigated BCJR algorithms whose calculation of recursions is limited to $M$ significant terms. The

8The math min phase model is obtained from 60 or more central samples of $h ( t )$ and begins with ≈ .031, .142, .344, . . .. To obtain a fair comparison, we have delayed the receiver in [11] by $K _ { p } = 1$ (see Section II) and set its first tap to 0. Its performance will otherwise be worse than shown in Fig. 13.   
9However, Douillard [4] reports that only 3 iterations are needed at 6–7 dB SNR.

application has been to simple ISI removal and turbo equalization of channels with spectral zeros and strong narrowband ISI, where $M$ is much smaller than the effective ISI state space. We have proposed several important innovations. An improvement to the minimum phase allpass filtering sharpens the focus of the ISI model energy. When combined with a delayed and slightly mismatched BCJR, the decoding allows a smaller $M$ without significant loss in BER. By adding a third lowcomplexity M-BCJR recursion, LLR quality is improved for practical values of $M$ , leading to a major BER improvement in turbo equalization. Other innovations are the use of single tentative soft symbol estimates to improve the reduced-trellis benchmark BCJR and a modified method for retaining backward recursion values. All these improvements work together to create a turbo equalizer of reasonable complexity, which in an FTN application can lead simultaneously to an energy saving of 4 dB and a bandwidth reduction of $3 5 \%$ compared to binary orthogonal signaling.

# ACKNOWLEDGMENTS

This work was supported by the Swedish Research Council (VR) through Grant 621-2003-3210, and by the Swedish Foundation for Strategic Research (SSF) through its Strategic Center for High Speed Wireless Communication at Lund.

# REFERENCES

[1] J. E. Mazo, “Faster-than-Nyquist signaling,” Bell Syst. Tech. J., vol. 54, pp. 1451–1462, Oct. 1975.   
[2] A. Liveris and C. N. Georghiades, “Exploiting faster-than-Nyquist signaling,” IEEE Trans. Commun., vol. 51, pp. 1502–1511, Sep. 2003.   
[3] F. Rusek and J. B. Anderson, “Multi-stream faster-than-Nyquist signaling,” IEEE Trans. Commun., vol. 57, pp. 1329–1340, May 2009.   
[4] C. Douillard et al., “Iterative correction of intersymbol interference: turbo equalization,” Eur. Trans. Telecommun., vol. 6, pp. 507–511, Sep./Oct. 1995.   
[5] K. Balachandran and J. B. Anderson, “Reduced complexity sequence detection for nonminimum phase intersymbol interference channels,” IEEE Trans. Inf. Theory, vol. 43, pp. 275–280, Jan. 1997.   
[6] V. Franz and J. B. Anderson, “Concatenated decoding with a reducedsearch BCJR algorithm,” IEEE J. Sel. Areas Commun., vol. 16, pp. 186– 195, Feb. 1998.   
[7] C. Fragouli, N. Seshadri, and W. Turin, “Reduced-trellis equalization using the BCJR algorithm,” Wireless Commun. & Mobile Computing, vol. 1, pp. 397–406, 2001.   
[8] G. Colavolpe, G. Ferrari, and R. Raheli, “Reduced-state BCJR type algorithms,” IEEE J. Sel. Areas Communs., vol. 19, pp. 848–859, May 2001.   
[9] M. Sikora and D. J. Costello, Jr., “A new SISO algorithm with application to turbo equalization,” in Proc. 2005 IEEE Int. Symp. Information Theory, pp. 2031–2035.   
[10] C. M. Vithanage, C. Andrieu, and R. J. Piechocki, “Novel reduced-state BCJR algorithms,” IEEE Trans. Commun., vol. 55, pp. 1144–1152, June 2007.   
[11] D. Fertonani, A. Barbieri, and G. Colavolpe, “Reduced-complexity BCJR algorithm for turbo equalization,” IEEE Trans. Commun., vol. 55, pp. 2279–2287, Dec. 2007.   
[12] D. Fertonani, A. Barbieri, and G. Colavolpe, “Novel graph-based algorithms for soft-output detection over dispersive channels,” in Proc. 2008 IEEE Global Communs. Conf.   
[13] R. R. Anderson and G. J. Foschini, “The minimum distance for MLSE digital data systems of limited complexity,” IEEE Trans. Inf. Theory, vol. 21, pp. 544–551, Sep. 1975.   
[14] A. Prlja, J. B. Anderson, and F. Rusek, “Receivers for faster-than-Nyquist signaling with and without turbo equalization,” in Proc. 2008 IEEE Int. Symp. Information Theory.

[15] J. B. Anderson, A. Prlja, and F. Rusek, “New reduced state space BCJR algorithms for the ISI channel,” in Proc. 2009 IEEE Int. Symp. Information Theory.   
[16] J. B. Anderson and M. Zeinali, “Best rate 1/2 convolutional codes for turbo equalization with severe ISI,” Proc. 2012 IEEE Int. Symp. Information Theory.   
[17] F. Rusek and J. B. Anderson, “Constrained capacities for faster-than-Nyquist signaling,” IEEE Trans. Inf. Theory, vol. 55, pp. 764–775, Feb. 2009.   
[18] J. B. Anderson, “Tree encoding of speech,” IEEE Trans. Inf. Theory, vol. 21, pp. 379–387, July 1975.   
[19] A. Duel-Hallen and C. Heegard, “Delayed decision-feedback sequence estimation,” IEEE Trans. Commun., vol. 37, pp. 428–436, May 1989.   
[20] M. Magarini, A. Spalvieri, and G. Tartara, “The mean-square delayed decision feedback sequence detector,” IEEE Trans. Commun., vol. 50, pp. 1462–1470, Sep. 2002.   
[21] U. L. Dang, W. H. Gerstacker, and S. T. M. Slock, “Maximum SINR prefiltering for reduced-state trellis-based equalization,” in Proc. 2011 IEEE Int. Conf. on Commun.   
[22] J. B. Anderson and A. Svensson, Coded Modulation Systems. Kluwer-Plenum, 2003.   
[23] S. ten Brink, “Convergence of iterative decoding,” IEE Electron. Lett., vol. 35, pp. 806–808, May 1999.   
[24] G. Colavolpe and A. Barbieri, “On MAP symbol detection for ISI channels using the Ungerboeck observation model,” IEEE Commun. Lett., vol. 9, pp. 720–722, Aug. 2005.

![](images/5e644b88d40d78e5617cd9fbf9bdf880c490ff6ce57d41492b3d1e84142b55a0.jpg)

Adnan Prlja was born in Banja Luka, Bosnia and Herzegovina in 1983. In March 2007 he received the degree of Master of Science (MSc) in Electrical Engineering from Lund University, Lund, Sweden, where he is currently finishing as a Ph.D. student at the Department of Electrical and Information Technology (EIT). His current research is mainly focused on iterative detection/decoding over channels with memory but it also includes other related topics from digital communications and information theory.

![](images/cfaf0b982523a90409469209b62b0a7a2f9babee067695f0155ee50992e24f49.jpg)

John B. Anderson was born in New York State in 1945. He received the B.S., M.S. and Ph.D. degrees in electrical engineering from Cornell University in 1967, 1969 and 1972. During 1972-80 he was on the electrical engineering faculty at McMaster University in Canada, and during 1981-98 he was Professor at Rensselaer Polytechnic Institute. Since 1998 he has held the Ericsson Chair in Digital Communication at Lund University, Sweden. He has held a number of visiting professorships in the US, Sweden, Canada and Germany. His research work

is in coding and communication algorithms, bandwidth-efficient coding, and applications of these to data transmission and compression. He has served widely as a consultant in these fields. Recently, he has been Director of the Swedish Strategic Research Foundation Center for High Speed Wireless Communication at Lund, and coordinator for Lund University of the Swedish national strategic funding scheme ELLIIT.

Dr. Anderson was a member of the IEEE Information Theory Society Board of Governors during 1980–87 and 2001–06, serving as the Society’s Vice-President and President (1985). In 1983 and 2006 he was Co-Chair of the IEEE International Symposium on Information Theory. In the IEEE publications sphere, he served on the Publications Board of IEEE on three occasions, and has been Editor-in-Chief of IEEE Press during 1994-96 and 2012–13. Since 1998 he has edited the Press book Series on Digital and Mobile Communication. He has also served as associate editor for several IEEE transactions.

Dr. Anderson is an author of six textbooks, including most recently Digital Transmission Engineering (IEEE Press, 2nd edition, 2005), Coded Modulation Systems (Plenum/Springer, 2003), and Understanding Information Transmission (IEEE Press, 2005). He is Fellow of the IEEE (1987) and received the Humboldt Research Prize (Germany) in 1991. In 1996 he was elected Swedish National Visiting Chair in Information Technology. He received the IEEE Third Millenium Medal in 2000.