use full_palette::ORANGE;
use image::DynamicImage;
use itertools::Itertools;
use plotters::{prelude::*, style::text_anchor::Pos};

use crate::datastructures::pwm::Pwm;

/// Nucleotides and their colors
const BASES: [&str; 4] = ["A", "C", "G", "T"];
const COLORS: [RGBColor; 4] = [RED, BLUE, GREEN, ORANGE];

const CACHE_DIR: &str = ".cache";

fn create_and_load_char_bitmaps() -> [DynamicImage; 4] {
    let mut images = Vec::new();
    for base in BASES {
        let file_path = format!("{}/{}.bmp", CACHE_DIR, base);
        if let Ok(true) = std::fs::exists(file_path.clone()) {
            images.push(image::open(&file_path).unwrap());
            continue;
        }

        if let Ok(false) = std::fs::exists(CACHE_DIR) {
            std::fs::create_dir(CACHE_DIR).unwrap();
        }
        let size_mul = 2;

        let root =
            BitMapBackend::new(&file_path, (75 * size_mul, 100 * size_mul)).into_drawing_area();
        root.fill(&WHITE).unwrap();

        let mut chart = ChartBuilder::on(&root)
            .margin(0)
            .build_cartesian_2d(0.0..1.0, 0.0..1.0)
            .unwrap();

        chart
            .draw_series(std::iter::once(Text::new(
                base,
                (0.5, 0.4),
                ("FiraCode Nerd Font Mono", 165 * size_mul)
                    .into_font()
                    .color(&COLORS[BASES.iter().position(|&x| x == base).unwrap()])
                    .pos(Pos::new(
                        plotters::style::text_anchor::HPos::Center,
                        plotters::style::text_anchor::VPos::Center,
                    )),
            )))
            .unwrap();

        chart
            .configure_mesh()
            .disable_mesh()
            .disable_axes()
            .draw()
            .unwrap();

        root.present().unwrap();

        images.push(image::open(&file_path).unwrap());
    }

    images.try_into().unwrap()
}

pub fn clear_cache() {
    if let Ok(true) = std::fs::exists(CACHE_DIR) {
        std::fs::remove_dir_all(CACHE_DIR).unwrap();
    }
}

pub fn plot_pwm(name: &str, pwm: &Pwm, score: f64) -> Result<(), Box<dyn std::error::Error>> {
    let images = create_and_load_char_bitmaps();
    let size_mul = 1;
    let root = BitMapBackend::new(name, (200 * pwm.len() as u32 * size_mul, 400 * size_mul))
        .into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .margin(10)
        .caption(format!("Sequence Logo: ({:.2})", score), ("sans-serif", 30))
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(0.0..pwm.len() as f64, 0.0..1.0)?; // Y-axis scaled for entropy

    chart
        .configure_mesh()
        .disable_mesh()
        .disable_axes()
        .draw()?;

    let (w, h) = chart.plotting_area().dim_in_pixel();

    for (i, column) in pwm.matrix.iter().enumerate() {
        let mut y_offset = 1.0;

        for (j, &freq) in column
            .iter()
            .enumerate()
            .sorted_by(|a, b| b.1.partial_cmp(a.1).unwrap())
        {
            let height = freq; // Scale frequencies
            let image = images[j].clone().resize_exact(
                w / pwm.len() as u32,
                (h as f64 * freq) as u32,
                image::imageops::FilterType::Gaussian,
            );

            let elem = BitMapElement::with_owned_buffer(
                (i as f64, y_offset),
                (w / pwm.len() as u32, (h as f64 * freq) as u32),
                image.into_bytes(),
            )
            .unwrap();

            chart.plotting_area().draw(&elem)?;

            y_offset -= height;
        }
    }

    root.present()?;
    println!("Sequence logo saved as {}", name);
    Ok(())
}
