use prodigal_rs::{predict, predict_meta, train, ProdigalConfig, Strand};

/// Read the sample FASTA and extract raw sequence bytes (stripping headers/newlines).
fn load_sample_sequences() -> Vec<(String, Vec<u8>)> {
    let fasta = std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("Prodigal/anthus_aco.fas"),
    )
    .unwrap();

    let mut seqs = Vec::new();
    let mut name = String::new();
    let mut seq = Vec::new();

    for line in fasta.lines() {
        if line.starts_with('>') {
            if !seq.is_empty() {
                seqs.push((name.clone(), seq.clone()));
                seq.clear();
            }
            name = line[1..].to_string();
        } else {
            seq.extend_from_slice(line.trim().as_bytes());
        }
    }
    if !seq.is_empty() {
        seqs.push((name, seq));
    }
    seqs
}

/// Concatenate all sequences into one (for training).
fn concat_sequences(seqs: &[(String, Vec<u8>)]) -> Vec<u8> {
    let mut all = Vec::new();
    for (_, s) in seqs {
        all.extend_from_slice(s);
    }
    all
}

#[test]
fn test_predict_meta_returns_genes() {
    let seqs = load_sample_sequences();
    // Pick a sequence that has genes
    for (name, seq) in &seqs {
        let genes = predict_meta(seq).unwrap();
        if genes.is_empty() {
            continue;
        }
        eprintln!(
            "[meta] {} ({} bp): {} genes predicted",
            name,
            seq.len(),
            genes.len()
        );

        // Verify basic properties
        for g in &genes {
            assert!(g.begin >= 1);
            assert!(g.end >= 1);
            assert!(g.begin != g.end);
            assert!(g.confidence >= 50.0);
            assert!(g.confidence <= 100.0);
            assert!(g.gc_content >= 0.0 && g.gc_content <= 1.0);
        }
        return;
    }
    panic!("No sequence produced genes");
}

#[test]
fn test_predict_meta_all_sequences() {
    let seqs = load_sample_sequences();
    let mut total_genes = 0;

    for (name, seq) in &seqs {
        let genes = predict_meta(seq).unwrap();
        total_genes += genes.len();
        eprintln!(
            "  {} ({} bp): {} genes",
            name, seq.len(), genes.len()
        );
    }
    eprintln!("[meta all] total: {} genes across {} sequences", total_genes, seqs.len());
    assert!(total_genes > 0);
}

#[test]
fn test_train_and_predict() {
    let seqs = load_sample_sequences();
    let all = concat_sequences(&seqs);

    // Train on the full dataset
    let training = train(&all).unwrap();
    eprintln!(
        "[train] GC={:.2}%, trans_table={}, uses_sd={}",
        training.gc() * 100.0,
        training.translation_table(),
        training.uses_sd()
    );

    // Predict on each sequence
    let mut total_genes = 0;
    for (name, seq) in &seqs {
        let genes = predict(seq, &training).unwrap();
        total_genes += genes.len();
        eprintln!(
            "  {} ({} bp): {} genes",
            name, seq.len(), genes.len()
        );
    }
    eprintln!("[single] total: {} genes", total_genes);
    assert!(total_genes > 0);
}

#[test]
fn test_training_save_load_roundtrip() {
    let seqs = load_sample_sequences();
    let all = concat_sequences(&seqs);

    let training = train(&all).unwrap();

    let dir = tempfile::TempDir::new().unwrap();
    let path = dir.path().join("test.trn");

    training.save(&path).unwrap();
    let loaded = prodigal_rs::TrainingData::load(&path).unwrap();

    assert_eq!(training.gc(), loaded.gc());
    assert_eq!(training.translation_table(), loaded.translation_table());
    assert_eq!(training.uses_sd(), loaded.uses_sd());

    // Predictions should be identical
    let seq = &seqs[0].1;
    let genes1 = predict(seq, &training).unwrap();
    let genes2 = predict(seq, &loaded).unwrap();

    assert_eq!(genes1.len(), genes2.len());
    for (g1, g2) in genes1.iter().zip(genes2.iter()) {
        assert_eq!(g1.begin, g2.begin);
        assert_eq!(g1.end, g2.end);
    }
}

#[test]
fn test_empty_sequence_error() {
    let result = predict_meta(b"");
    assert!(result.is_err());
}

#[test]
fn test_short_sequence_training_error() {
    let result = train(b"ATGCGATCG");
    assert!(result.is_err());
}

#[test]
fn test_gene_properties() {
    let seqs = load_sample_sequences();
    for (_, seq) in &seqs {
        let genes = predict_meta(seq).unwrap();
        if genes.is_empty() {
            continue;
        }
        let g = &genes[0];

        // Check strand
        assert!(g.strand == Strand::Forward || g.strand == Strand::Reverse);

        // Check RBS motif is non-empty
        assert!(!g.rbs_motif.is_empty());

        // Check scores are finite
        assert!(g.score.is_finite());
        assert!(g.cscore.is_finite());
        assert!(g.sscore.is_finite());
        return;
    }
}

#[test]
fn test_custom_config() {
    let seqs = load_sample_sequences();
    let seq = &seqs[0].1; // use first sequence

    let config = ProdigalConfig {
        closed_ends: true,
        ..Default::default()
    };
    let genes = prodigal_rs::predict_meta_with(seq, &config).unwrap();
    // With closed ends, edge genes are suppressed
    for g in &genes {
        assert!(!g.partial.0 && !g.partial.1, "closed_ends should prevent partial genes");
    }
}
