use std::fs;
use std::error::Error;

fn initialize_score_matrix(m: usize, n: usize, gap_penalty: i32) -> Vec<Vec<i32>> {
    let mut score = vec![vec![0; n + 1]; m + 1];

    // Erste Spalte (Lücken in seq2)
    for i in 0..=m {
        score[i][0] = (i as i32) * gap_penalty;
    }

    // Erste Zeile (Lücken in seq1)
    for j in 0..=n {
        score[0][j] = (j as i32) * gap_penalty;
    }

    score
}

fn calculate_match_score(char1: char, char2: char, match_score: i32, mismatch_penalty: i32)
    -> i32 {
    if char1 == char2 {
        match_score
    } else {
        mismatch_penalty
    }
}

fn fill_score_matrix(
    score: &mut Vec<Vec<i32>>,
    s1: &[char],
    s2: &[char],
    match_score: i32,
    mismatch_penalty: i32,
    gap_penalty: i32,
) {
    let m = s1.len();
    let n = s2.len();

    for i in 1..=m {
        for j in 1..=n {
            let match_mismatch = calculate_match_score(s1[i - 1], s2[j - 1],
                                                       match_score, mismatch_penalty);

            score[i][j] = [
                score[i - 1][j - 1] + match_mismatch, // Diagonal
                score[i - 1][j] + gap_penalty,        // Oben (Lücke in seq2)
                score[i][j - 1] + gap_penalty,        // Links (Lücke in seq1)
            ]
                .iter()
                .max()
                .copied()
                .unwrap();
        }
    }
}

fn traceback_alignment(
    score: &Vec<Vec<i32>>,
    s1: &[char],
    s2: &[char],
    match_score: i32,
    mismatch_penalty: i32,
    gap_penalty: i32,
) -> (String, String) {
    let mut aligned1 = String::new();
    let mut aligned2 = String::new();
    let (mut i, mut j) = (s1.len(), s2.len());

    while i > 0 || j > 0 {
        if i > 0 && j > 0 && score[i][j] == score[i - 1][j - 1] +
            calculate_match_score(s1[i - 1], s2[j - 1], match_score,
                                  mismatch_penalty) {
            // Diagonal: Match oder Mismatch
            aligned1.insert(0, s1[i - 1]);
            aligned2.insert(0, s2[j - 1]);
            i -= 1;
            j -= 1;
        } else if i > 0 && score[i][j] == score[i - 1][j] + gap_penalty {
            // Oben: Lücke in seq2
            aligned1.insert(0, s1[i - 1]);
            aligned2.insert(0, '-');
            i -= 1;
        } else {
            // Links: Lücke in seq1
            aligned1.insert(0, '-');
            aligned2.insert(0, s2[j - 1]);
            j -= 1;
        }
    }

    (aligned1, aligned2)
}

fn print_alignment(aligned1: &str, aligned2: &str, score: i32) {
    let mut match_line = String::new();

    for (c1, c2) in aligned1.chars().zip(aligned2.chars()) {
        if c1 == c2 {
            match_line.push('|');
        } else if c1 == '-' || c2 == '-' {
            match_line.push(' ');
        } else {
            match_line.push('x');
        }
    }

    println!("Sequence 1: {}", aligned1);
    println!("            {}", match_line);
    println!("Sequence 2: {}", aligned2);
    println!("Score:      {}", score);
}

fn needleman_wunsch(seq1: &str, seq2: &str, match_score: i32, mismatch_penalty: i32,
                    gap_penalty: i32) -> (String, String, i32) {
    let m = seq1.len();
    let n = seq2.len();
    let s1: Vec<char> = seq1.chars().collect();
    let s2: Vec<char> = seq2.chars().collect();

    // Score-Matrix initialisieren
    let mut score = initialize_score_matrix(m, n, gap_penalty);

    // Score-Matrix auffüllen
    fill_score_matrix(&mut score, &s1, &s2, match_score, mismatch_penalty, gap_penalty);

    // Traceback durchführen
    let (aligned1, aligned2) = traceback_alignment(&score, &s1, &s2, match_score,
                                                   mismatch_penalty, gap_penalty);

    (aligned1, aligned2, score[m][n])
}

fn main() -> Result<(), Box<dyn Error>> {
    let seq1 = fs::read_to_string("src/data/HBB_wildtype.txt")?;
    let seq2 = fs::read_to_string("src/data/HBB_sickle_cell_allele.txt")?;
    let (aligned1, aligned2, score) = needleman_wunsch(&seq1, &seq2,
                                                       1, -1, -1);

    print_alignment(&aligned1, &aligned2, score);
    Ok(())
}
