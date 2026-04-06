use crate::types::{MetagenomicBin, Training};
use std::os::raw::c_int;

extern "C" {
    fn initialize_metagenome_0(tptr: *mut Training);
    fn initialize_metagenome_1(tptr: *mut Training);
    fn initialize_metagenome_2(tptr: *mut Training);
    fn initialize_metagenome_3(tptr: *mut Training);
    fn initialize_metagenome_4(tptr: *mut Training);
    fn initialize_metagenome_5(tptr: *mut Training);
    fn initialize_metagenome_6(tptr: *mut Training);
    fn initialize_metagenome_7(tptr: *mut Training);
    fn initialize_metagenome_8(tptr: *mut Training);
    fn initialize_metagenome_9(tptr: *mut Training);
    fn initialize_metagenome_10(tptr: *mut Training);
    fn initialize_metagenome_11(tptr: *mut Training);
    fn initialize_metagenome_12(tptr: *mut Training);
    fn initialize_metagenome_13(tptr: *mut Training);
    fn initialize_metagenome_14(tptr: *mut Training);
    fn initialize_metagenome_15(tptr: *mut Training);
    fn initialize_metagenome_16(tptr: *mut Training);
    fn initialize_metagenome_17(tptr: *mut Training);
    fn initialize_metagenome_18(tptr: *mut Training);
    fn initialize_metagenome_19(tptr: *mut Training);
    fn initialize_metagenome_20(tptr: *mut Training);
    fn initialize_metagenome_21(tptr: *mut Training);
    fn initialize_metagenome_22(tptr: *mut Training);
    fn initialize_metagenome_23(tptr: *mut Training);
    fn initialize_metagenome_24(tptr: *mut Training);
    fn initialize_metagenome_25(tptr: *mut Training);
    fn initialize_metagenome_26(tptr: *mut Training);
    fn initialize_metagenome_27(tptr: *mut Training);
    fn initialize_metagenome_28(tptr: *mut Training);
    fn initialize_metagenome_29(tptr: *mut Training);
    fn initialize_metagenome_30(tptr: *mut Training);
    fn initialize_metagenome_31(tptr: *mut Training);
    fn initialize_metagenome_32(tptr: *mut Training);
    fn initialize_metagenome_33(tptr: *mut Training);
    fn initialize_metagenome_34(tptr: *mut Training);
    fn initialize_metagenome_35(tptr: *mut Training);
    fn initialize_metagenome_36(tptr: *mut Training);
    fn initialize_metagenome_37(tptr: *mut Training);
    fn initialize_metagenome_38(tptr: *mut Training);
    fn initialize_metagenome_39(tptr: *mut Training);
    fn initialize_metagenome_40(tptr: *mut Training);
    fn initialize_metagenome_41(tptr: *mut Training);
    fn initialize_metagenome_42(tptr: *mut Training);
    fn initialize_metagenome_43(tptr: *mut Training);
    fn initialize_metagenome_44(tptr: *mut Training);
    fn initialize_metagenome_45(tptr: *mut Training);
    fn initialize_metagenome_46(tptr: *mut Training);
    fn initialize_metagenome_47(tptr: *mut Training);
    fn initialize_metagenome_48(tptr: *mut Training);
    fn initialize_metagenome_49(tptr: *mut Training);
}

static DESCS: [(&[u8], &[u8], f64); 50] = [
    (b"Mycoplasma_bovis_PG45\0", b"B\0", 29.31),
    (b"Mycoplasma_pneumoniae_M129\0", b"B\0", 40.01),
    (b"Mycoplasma_suis_Illinois\0", b"B\0", 31.08),
    (b"Aeropyrum_pernix_K1\0", b"A\0", 56.31),
    (b"Akkermansia_muciniphila_ATCC_BAA_835\0", b"B\0", 55.76),
    (b"Anaplasma_marginale_Maries\0", b"B\0", 49.76),
    (b"Anaplasma_phagocytophilum_HZ\0", b"B\0", 41.64),
    (b"Archaeoglobus_fulgidus_DSM_4304\0", b"A\0", 48.58),
    (b"Bacteroides_fragilis_NCTC_9343\0", b"B\0", 43.19),
    (b"Brucella_canis_ATCC_23365\0", b"B\0", 57.21),
    (b"Burkholderia_rhizoxinica_HKI_454\0", b"B\0", 59.70),
    (b"Candidatus_Amoebophilus_asiaticus_5a2\0", b"B\0", 35.05),
    (b"Candidatus_Korarchaeum_cryptofilum_OPF8\0", b"A\0", 49.00),
    (b"Catenulispora_acidiphila_DSM_44928\0", b"B\0", 69.77),
    (b"Cenarchaeum_symbiosum_B\0", b"A\0", 57.19),
    (b"Chlorobium_phaeobacteroides_BS1\0", b"B\0", 48.93),
    (b"Chlorobium_tepidum_TLS\0", b"B\0", 56.53),
    (b"Desulfotomaculum_acetoxidans_DSM_771\0", b"B\0", 41.55),
    (b"Desulfurococcus_kamchatkensis_1221n\0", b"B\0", 45.34),
    (b"Erythrobacter_litoralis_HTCC2594\0", b"B\0", 63.07),
    (b"Escherichia_coli_UMN026\0", b"B\0", 50.72),
    (b"Haloquadratum_walsbyi_DSM_16790\0", b"A\0", 47.86),
    (b"Halorubrum_lacusprofundi_ATCC_49239\0", b"A\0", 57.14),
    (b"Hyperthermus_butylicus_DSM_5456\0", b"A\0", 53.74),
    (b"Ignisphaera_aggregans_DSM_17230\0", b"A\0", 35.69),
    (b"Marinobacter_aquaeolei_VT8\0", b"B\0", 57.27),
    (b"Methanopyrus_kandleri_AV19\0", b"A\0", 61.16),
    (b"Methanosphaerula_palustris_E1_9c\0", b"A\0", 55.35),
    (b"Methanothermobacter_thermautotrophicus_Delta_H\0", b"B\0", 49.54),
    (b"Methylacidiphilum_infernorum_V4\0", b"B\0", 45.48),
    (b"Mycobacterium_leprae_TN\0", b"B\0", 57.80),
    (b"Natrialba_magadii_ATCC_43099\0", b"A\0", 61.42),
    (b"Orientia_tsutsugamushi_Boryong\0", b"B\0", 30.53),
    (b"Pelotomaculum_thermopropionicum_SI\0", b"B\0", 52.96),
    (b"Prochlorococcus_marinus_MIT_9313\0", b"B\0", 50.74),
    (b"Pyrobaculum_aerophilum_IM2\0", b"A\0", 51.36),
    (b"Ralstonia_solanacearum_PSI07\0", b"B\0", 66.13),
    (b"Rhizobium_NGR234\0", b"B\0", 58.49),
    (b"Rhodococcus_jostii_RHA1\0", b"B\0", 65.05),
    (b"Rickettsia_conorii_Malish_7\0", b"B\0", 32.44),
    (b"Rothia_dentocariosa_ATCC_17931\0", b"B\0", 53.69),
    (b"Shigella_dysenteriae_Sd197\0", b"B\0", 51.25),
    (b"Synechococcus_CC9605\0", b"B\0", 59.22),
    (b"Synechococcus_JA_2_3B_a_2_13_\0", b"B\0", 58.45),
    (b"Thermoplasma_volcanium_GSS1\0", b"A\0", 39.92),
    (b"Treponema_pallidum_Nichols\0", b"B\0", 52.77),
    (b"Tropheryma_whipplei_TW08_27\0", b"B\0", 46.31),
    (b"Xenorhabdus_nematophila_ATCC_19061\0", b"B\0", 44.15),
    (b"Xylella_fastidiosa_Temecula1\0", b"B\0", 51.78),
    (b"_Nostoc_azollae__0708\0", b"B\0", 38.45),
];

static INIT_FNS: [unsafe extern "C" fn(*mut Training); 50] = [
    initialize_metagenome_0,  initialize_metagenome_1,
    initialize_metagenome_2,  initialize_metagenome_3,
    initialize_metagenome_4,  initialize_metagenome_5,
    initialize_metagenome_6,  initialize_metagenome_7,
    initialize_metagenome_8,  initialize_metagenome_9,
    initialize_metagenome_10, initialize_metagenome_11,
    initialize_metagenome_12, initialize_metagenome_13,
    initialize_metagenome_14, initialize_metagenome_15,
    initialize_metagenome_16, initialize_metagenome_17,
    initialize_metagenome_18, initialize_metagenome_19,
    initialize_metagenome_20, initialize_metagenome_21,
    initialize_metagenome_22, initialize_metagenome_23,
    initialize_metagenome_24, initialize_metagenome_25,
    initialize_metagenome_26, initialize_metagenome_27,
    initialize_metagenome_28, initialize_metagenome_29,
    initialize_metagenome_30, initialize_metagenome_31,
    initialize_metagenome_32, initialize_metagenome_33,
    initialize_metagenome_34, initialize_metagenome_35,
    initialize_metagenome_36, initialize_metagenome_37,
    initialize_metagenome_38, initialize_metagenome_39,
    initialize_metagenome_40, initialize_metagenome_41,
    initialize_metagenome_42, initialize_metagenome_43,
    initialize_metagenome_44, initialize_metagenome_45,
    initialize_metagenome_46, initialize_metagenome_47,
    initialize_metagenome_48, initialize_metagenome_49,
];

#[no_mangle]
pub unsafe extern "C" fn initialize_metagenomic_bins(meta: *mut MetagenomicBin) {
    for i in 0..50 {
        let m = &mut *meta.add(i);
        INIT_FNS[i](m.tinf);
        libc::sprintf(
            m.desc.as_mut_ptr(),
            b"%d|%s|%s|%.1f|%d|%d\0".as_ptr() as *const libc::c_char,
            i as c_int,
            DESCS[i].0.as_ptr() as *const libc::c_char,
            DESCS[i].1.as_ptr() as *const libc::c_char,
            DESCS[i].2,
            (*m.tinf).trans_table,
            (*m.tinf).uses_sd,
        );
    }
}
