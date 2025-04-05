use clap::Parser;
use std::path::Path;

mod args;
mod datastructures;
mod fasta_reader;
mod motif_finder;
mod plot;

fn main() {
    let args = args::Args::parse();
    let path = Path::new(&args.input_file);

    if !path.exists() {
        panic!("File does not exist: {:?}", path);
    }

    let seqs = fasta_reader::read_fasta(path);

    println!("Input file: {:?}", path);
    println!("Sequences: {}", seqs.len());
    println!(
        "Average sequence length: {:.2}",
        seqs.iter().map(|x| x.len()).sum::<usize>() as f64 / seqs.len() as f64
    );
    print!("====================");

    motif_finder::motif_finder(seqs, args.plot_sequence_logos, args.hits);
}
