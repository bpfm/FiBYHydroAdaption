pro mass_fractions

POPIII_N_MASS_BIN = 200

POPIII_IMF_Exponent = 2.35

POPIII_IMF_MinMass_MSUN =  21.0 ; solar masses
POPIII_IMF_MaxMass_MSUN = 500.0 ; solar masses

; allocate arrays
popiii_imf_mass_bin = dblarr(POPIII_N_MASS_BIN)
popiii_imf_mass_bin_log10 = dblarr(POPIII_N_MASS_BIN)
popiii_imf_by_number = dblarr(POPIII_N_MASS_BIN)

; define bins and IMF
lm_min = alog10(POPIII_IMF_MinMass_MSUN)
lm_max = alog10(POPIII_IMF_MaxMass_MSUN)

dlm = (lm_max - lm_min) / (POPIII_N_MASS_BIN - 1)

for i = 0, POPIII_N_MASS_BIN - 1 do begin
    lmass = lm_min + i * dlm    
    mass = 10.0^lmass

    popiii_imf_by_number[i] = 1.0 / mass^POPIII_IMF_Exponent
    popiii_imf_mass_bin[i] = mass
    popiii_imf_mass_bin_log10[i] = lmass
endfor

; normalize integral of IMF to 1
normalization = total(popiii_imf_by_number * popiii_imf_mass_bin * popiii_imf_mass_bin)

normalization -= 0.5 * popiii_imf_by_number[0] * $
  popiii_imf_mass_bin[0] * popiii_imf_mass_bin[0]

normalization -= 0.5 * popiii_imf_by_number[POPIII_N_MASS_BIN-1] * $
  popiii_imf_mass_bin[POPIII_N_MASS_BIN-1] * popiii_imf_mass_bin[POPIII_N_MASS_BIN-1]

normalization *= dlm * alog(10.0)

print, 'Total IMF integral = ', normalization

popiii_imf_by_number /= normalization

print, 'Total IMF integral = ', (total(popiii_imf_by_number * $
                                       popiii_imf_mass_bin * $
                                       popiii_imf_mass_bin) - $
                                 0.5 * popiii_imf_by_number[0] * popiii_imf_mass_bin[0] * popiii_imf_mass_bin[0] - $
                                 0.5 * popiii_imf_by_number[POPIII_N_MASS_BIN-1] * popiii_imf_mass_bin[POPIII_N_MASS_BIN-1] * popiii_imf_mass_bin[POPIII_N_MASS_BIN-1] ) * $
  dlm * alog(10.0)

; ; test Hydrogen fraction
; yield = rdtab2('/data/mpe/caius/POPIII/TEST_YIELDS/POPIIIyields.Hydrogen.spline',cols=[1,2])

; yield_bin = yield[1,*]

; mas_fraction = (total(popiii_imf_by_number * $
;                       popiii_imf_mass_bin * $
;                       yield_bin) - $
;                 0.5 * popiii_imf_by_number[0] * popiii_imf_mass_bin[0] * yield_bin[0] - $
;                 0.5 * popiii_imf_by_number[POPIII_N_MASS_BIN-1] * popiii_imf_mass_bin[POPIII_N_MASS_BIN-1] * yield_bin[POPIII_N_MASS_BIN-1] ) * $
;   dlm * alog(10.0)

; print, 'H mass fraction = ', mas_fraction

; stop

; ; test Carbon fraction
; yield = rdtab2('/data/mpe/caius/POPIII/TEST_YIELDS/POPIIIyields.Carbon.spline',cols=[1,2])

; yield_bin = yield[1,*]

; mas_fraction = (total(popiii_imf_by_number * $
;                       popiii_imf_mass_bin * $
;                       yield_bin) - $
;                 0.5 * popiii_imf_by_number[0] * popiii_imf_mass_bin[0] * yield_bin[0] - $
;                 0.5 * popiii_imf_by_number[POPIII_N_MASS_BIN-1] * popiii_imf_mass_bin[POPIII_N_MASS_BIN-1] * yield_bin[POPIII_N_MASS_BIN-1] ) * $
;   dlm * alog(10.0)

; print, 'C mass fraction = ', mas_fraction

; stop

; ; test Nitrogen fraction
; yield = rdtab2('/data/mpe/caius/POPIII/TEST_YIELDS/POPIIIyields.Nitrogen.spline',cols=[1,2])

; yield_bin = yield[1,*]

; mas_fraction = (total(popiii_imf_by_number * $
;                       popiii_imf_mass_bin * $
;                       yield_bin) - $
;                 0.5 * popiii_imf_by_number[0] * popiii_imf_mass_bin[0] * yield_bin[0] - $
;                 0.5 * popiii_imf_by_number[POPIII_N_MASS_BIN-1] * popiii_imf_mass_bin[POPIII_N_MASS_BIN-1] * yield_bin[POPIII_N_MASS_BIN-1] ) * $
;   dlm * alog(10.0)

; print, 'N mass fraction = ', mas_fraction

; stop

; ; test Oxygen fraction
; yield = rdtab2('/data/mpe/caius/POPIII/TEST_YIELDS/POPIIIyields.Oxygen.spline',cols=[1,2])

; yield_bin = yield[1,*]

; mas_fraction = (total(popiii_imf_by_number * $
;                       popiii_imf_mass_bin * $
;                       yield_bin) - $
;                 0.5 * popiii_imf_by_number[0] * popiii_imf_mass_bin[0] * yield_bin[0] - $
;                 0.5 * popiii_imf_by_number[POPIII_N_MASS_BIN-1] * popiii_imf_mass_bin[POPIII_N_MASS_BIN-1] * yield_bin[POPIII_N_MASS_BIN-1] ) * $
;   dlm * alog(10.0)

; print, 'O mass fraction = ', mas_fraction

; stop

; black hole mass fraction [260-500] solar masses
ejecta = rdtab2('/data/mpe/caius/POPIII/TEST_YIELDS/POPIIIejecta.spline',cols=[1,2])

ejected_mass = (total(popiii_imf_by_number * $
                      popiii_imf_mass_bin * $
                      ejecta[1,*]) - $
                0.5 * popiii_imf_by_number[0] * popiii_imf_mass_bin[0] * ejecta[1,0] - $
                0.5 * popiii_imf_by_number[POPIII_N_MASS_BIN-1] * popiii_imf_mass_bin[POPIII_N_MASS_BIN-1] * ejecta[1,POPIII_N_MASS_BIN-1] ) * $
  dlm * alog(10.0)

total_mass = (total(popiii_imf_by_number * $
                    popiii_imf_mass_bin * $
                    popiii_imf_mass_bin) - $
              0.5 * popiii_imf_by_number[0] * popiii_imf_mass_bin[0] * popiii_imf_mass_bin[0] - $
              0.5 * popiii_imf_by_number[POPIII_N_MASS_BIN-1] * popiii_imf_mass_bin[POPIII_N_MASS_BIN-1] * popiii_imf_mass_bin[POPIII_N_MASS_BIN-1] ) * $
  dlm * alog(10.0)

print, 'Total ejected mass = ', ejected_mass, '   total mass = ', total_mass

mass_min = 26
mass_max = 100

index = where(popiii_imf_mass_bin le mass_min)
index_min = index[n_elements(index)-1]

index = where(popiii_imf_mass_bin ge mass_max)
index_max = index[0]

print,index_min,index_max

ejected_mass_25 = (total(popiii_imf_by_number[index_min:index_max] * $
                         popiii_imf_mass_bin[index_min:index_max] * $
                         ejecta[1,index_min:index_max]) - $
                0.5 * popiii_imf_by_number[index_min] * popiii_imf_mass_bin[index_min] * ejecta[1,index_min] - $
                0.5 * popiii_imf_by_number[index_max] * popiii_imf_mass_bin[index_max] * ejecta[1,index_max] ) * $
  dlm * alog(10.0)

total_mass_25 = (total(popiii_imf_by_number[index_min:index_max] * $
                       popiii_imf_mass_bin[index_min:index_max] * $
                       popiii_imf_mass_bin[index_min:index_max]) - $
                 0.5 * popiii_imf_by_number[index_min] * popiii_imf_mass_bin[index_min] * popiii_imf_mass_bin[index_min] - $
                 0.5 * popiii_imf_by_number[index_max] * popiii_imf_mass_bin[index_max] * popiii_imf_mass_bin[index_max] ) * $
  dlm * alog(10.0)

print, 'BH mass fraction = ', total_mass_25 - ejected_mass_25, total_mass_25, ejected_mass_25

mass_min = 260
mass_max = 500

index = where(popiii_imf_mass_bin le mass_min)
index_min = index[n_elements(index)-1]

index = where(popiii_imf_mass_bin ge mass_max)
index_max = index[0]

popiii_bh_fraction = (total(popiii_imf_by_number[index_min:index_max] * $
                            popiii_imf_mass_bin[index_min:index_max] * $
                            popiii_imf_mass_bin[index_min:index_max]) - $
                      0.5 * popiii_imf_by_number[index_min] * popiii_imf_mass_bin[index_min] * popiii_imf_mass_bin[index_min] - $
                      0.5 * popiii_imf_by_number[index_max] * popiii_imf_mass_bin[index_max] * popiii_imf_mass_bin[index_max] ) * $
  dlm * alog(10.0)

print, 'BH mass fraction = ', popiii_bh_fraction






































end
