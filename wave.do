onerror {resume}
quietly WaveActivateNextPane {} 0
add wave -noupdate /FFTtest/FFT/A
add wave -noupdate /FFTtest/FFT/B
update
WaveRestoreZoom {0 ps} {74100 ps}
run -all