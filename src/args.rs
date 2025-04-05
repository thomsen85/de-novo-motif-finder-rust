use clap::Parser;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Args {
    /// Path of the input file
    #[arg(short, long)]
    pub input_file: String,

    /// If you want to create sequence logos
    #[arg(short, long)]
    pub plot_sequence_logos: bool,

    /// Amount of motifs to return
    #[arg(long, default_value = "3")]
    pub hits: usize,

    /// Max length of the motifs
    #[arg(long, default_value = "20")]
    pub max_length: usize,

    /// Min length of the motifs
    #[arg(long, default_value = "6")]
    pub min_length: usize,

    /// Allowed gap size in motifs
    #[arg(long, default_value = "1")]
    pub gap_size: usize,
}
