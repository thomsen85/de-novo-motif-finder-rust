use std::path::Path;

use clap::Parser;

mod bio;
mod datastructures;
mod fasta_reader;
mod motif_finder;
mod plot;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Name of the person to greet
    #[arg(long)]
    input_file: String,

    /// Amount of hits to return
    #[arg(long, default_value = "3")]
    hits: usize,

    /// If you want to create sequence logos
    #[arg(long)]
    create_logos: bool,
}

fn main() {
    let args = Args::parse();
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

    motif_finder::motif_finder(seqs, args.create_logos, args.hits);
}
