ipp_1_csv = "data/ipp-1-threads.csv"
ipp_2_csv = "data/ipp-2-threads.csv"
ipp_4_csv = "data/ipp-4-threads.csv"
ipp_8_csv = "data/ipp-8-threads.csv"

##############################################################
# common options
#
# recall that
# lc - linecolor
# lt - linetype
# lw - linewidth
# pt - pointtype
# ps - pointsize
# w  - with
# lp - linespoints
# ls - linestyle
##############################################################

set terminal pdfcairo font "Roboto Sans,12" size 4,3 background rgb 'white'
set termoption enhanced # for superscripts
set datafile separator ","

set xtics 2 nomirror
set ytics nomirror

set xrange [*:*]

set grid back lt 1 dt 3 lc rgb 'grey'

set border 3 back

set style line 1 lc rgb "#000000" linewidth 2.5 pointtype 5 pointsize 0.85
set style line 2 lc rgb "#D53E4F" linewidth 2.5 pointtype 9 pointsize 0.85
set style line 4 lc rgb "#FFA500" linewidth 2.5 pointtype 7 pointsize 0.85
set style line 8 lc rgb "#008000" linewidth 2.5 pointtype 11 pointsize 0.85


set key top left
set key spacing 1.5
set key samplen 4

set logscale x 2
set logscale y 10
set format x "2^{%L}"; # disregards mantissa
set format y "10^{%L}"; # disregards mantissa

##############################################################
##############################################################
##############################################################

set output "plots/ipp-1-thread.pdf"
set xlabel "Number of inputs" offset 3,0,0
set ylabel "Time (s)" offset -1,0,0

plot \
     ipp_1_csv using (column('size')):(column('direct')) title 'Direct' with lp ls 1, \
     ipp_1_csv using (column('size')):(column('prover')) title 'Prover' with lp ls 2, \
     ipp_1_csv using (column('size')):(column('verifier')) title 'Verifier' with lp ls 4, \

##############################################################
##############################################################

set output "plots/ipp-direct-across-threads.pdf"
set xlabel "Number of inputs" offset 3,0,0
set ylabel "Time to directly compute pairing product (s)" offset -1,0,0

plot \
     ipp_1_csv using (column('size')):(column('direct')) title '1 thread' with lp ls 1, \
     ipp_2_csv using (column('size')):(column('direct')) title '2 threads' with lp ls 2, \
     ipp_4_csv using (column('size')):(column('direct')) title '4 threads' with lp ls 4, \
     ipp_8_csv using (column('size')):(column('direct')) title '8 threads' with lp ls 8, \

##############################################################
##############################################################

set output "plots/ipp-prover-across-threads.pdf"
set xlabel "Number of inputs" offset 3,0,0
set ylabel "Prover time (s)" offset -1,0,0

plot \
     ipp_1_csv using (column('size')):(column('prover')) title '1 thread' with lp ls 1, \
     ipp_2_csv using (column('size')):(column('prover')) title '2 threads' with lp ls 2, \
     ipp_4_csv using (column('size')):(column('prover')) title '4 threads' with lp ls 4, \
     ipp_8_csv using (column('size')):(column('prover')) title '8 threads' with lp ls 8, \

##############################################################
##############################################################

set output "plots/ipp-verifier-across-threads.pdf"
set xlabel "Number of inputs" offset 3,0,0
set ylabel "Verifier time (s)" offset -1,0,0

plot \
     ipp_1_csv using (column('size')):(column('verifier')) title '1 thread' with lp ls 1, \
     ipp_2_csv using (column('size')):(column('verifier')) title '2 threads' with lp ls 2, \
     ipp_4_csv using (column('size')):(column('verifier')) title '4 threads' with lp ls 4, \
     ipp_8_csv using (column('size')):(column('verifier')) title '8 threads' with lp ls 8, \

##############################################################
##############################################################

unset logscale y
unset format y
set output "plots/ipp-ratio.pdf"
set xlabel "Number of inputs" offset 3,0,0
set ylabel "Direct product time vs verifier time" offset -1,0,0

plot \
     ipp_1_csv using (column('size')):(column('direct')/column('verifier')) title '1 thread' with lp ls 1, \
     ipp_2_csv using (column('size')):(column('direct')/column('verifier')) title '2 threads' with lp ls 2, \
     ipp_4_csv using (column('size')):(column('direct')/column('verifier')) title '4 threads' with lp ls 4, \
     ipp_8_csv using (column('size')):(column('direct')/column('verifier')) title '8 threads' with lp ls 8, \

##############################################################
##############################################################

set output "plots/ipp-1-thread-linear.pdf"
set xlabel "Number of inputs" offset 3,0,0
set ylabel "Time (s)" offset -1,0,0

unset logscale x
set xtics (2, 2**14, 2**15, 2**16, 2**17)
plot \
     "<(sed '3,14d' data/ipp-1-threads.csv)" using (column('size')):(column('direct')) title 'Direct' with lp ls 1, \
     "<(sed '3,14d' data/ipp-1-threads.csv)" using (column('size')):(column('prover')) title 'Prover' with lp ls 2, \
     "<(sed '3,14d' data/ipp-1-threads.csv)" using (column('size')):(column('verifier')) title 'Verifier' with lp ls 4, \
