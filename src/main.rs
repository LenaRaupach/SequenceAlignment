fn needleman_wunsch(seq1: &str, seq2: &str, match_score: i32, mismatch_penalty: i32, gap_penalty: i32) -> (String, String, i32) {
    let m = seq1.len();
    let n = seq2.len();
    let s1: Vec<char> = seq1.chars().collect();
    let s2: Vec<char> = seq2.chars().collect();

    // Dynamische Programmierung: Score-Matrix initialisieren
    let mut score = vec![vec![0; n + 1]; m + 1];

    for i in 0..=m {
        score[i][0] = (i as i32) * gap_penalty;
    }
    for j in 0..=n {
        score[0][j] = (j as i32) * gap_penalty;
    }

    // Auffüllen der Score-Matrix
    for i in 1..=m {
        for j in 1..=n {
            let match_mismatch = if s1[i - 1] == s2[j - 1] {
                match_score
            } else {
                mismatch_penalty
            };
            score[i][j] = *[
                score[i - 1][j - 1] + match_mismatch, // Diagonal
                score[i - 1][j] + gap_penalty,        // Oben (Lücke in seq2)
                score[i][j - 1] + gap_penalty,        // Links (Lücke in seq1)
            ]
                .iter()
                .max()
                .unwrap();
        }
    }

    // Traceback: Alignment rekonstruieren
    let mut aligned1 = String::new();
    let mut aligned2 = String::new();
    let (mut i, mut j) = (m, n);

    while i > 0 || j > 0 {
        if i > 0 && j > 0 && score[i][j] == score[i - 1][j - 1] +
            if s1[i - 1] == s2[j - 1] { match_score } else { mismatch_penalty } {
            aligned1.insert(0, s1[i - 1]);
            aligned2.insert(0, s2[j - 1]);
            i -= 1;
            j -= 1;
        } else if i > 0 && score[i][j] == score[i - 1][j] + gap_penalty {
            aligned1.insert(0, s1[i - 1]);
            aligned2.insert(0, '-');
            i -= 1;
        } else {
            aligned1.insert(0, '-');
            aligned2.insert(0, s2[j - 1]);
            j -= 1;
        }
    }

    (aligned1, aligned2, score[m][n])
}

fn main() {
    let seq1 = "GATTACA";
    let seq2 = "GCATGCU";
    let (aligned1, aligned2, score) = needleman_wunsch(seq1, seq2, 1, -1, -1);

    println!("Alignment 1: {}", aligned1);
    println!("Alignment 2: {}", aligned2);
    println!("Score: {}", score);
}
