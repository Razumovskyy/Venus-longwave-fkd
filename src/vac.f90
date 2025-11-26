    !   - Uses: kd_data_mod, planck_mod, co2_corr_mod
    !   - Exposes: compute_vac_for_atm(za, wa, ta, gas_names, ro_co2, ro_h2o, ro_so2, vac_out, pl_out)
    !   - Calls into kd loaders once per run, applies CO₂ T-correction, interpolates KD to levels, writes VAC/
    !     PL arrays.