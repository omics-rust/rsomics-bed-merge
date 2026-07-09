use std::io;
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_bed_merge::{merge, merge_stdin};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-bed-merge", disable_help_flag = true)]
pub struct Cli {
    /// Input sorted BED (default: stdin)
    input: Option<PathBuf>,
    /// Output BED (default: stdout)
    #[arg(short = 'o', long)]
    output: Option<PathBuf>,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let mut stdout_lock;
        let mut file_out;
        let mut sink;
        // Under --json the framework owns stdout with the result envelope, so the
        // BED stream must not be written there too. Route it to -o when given,
        // otherwise to a sink so the envelope is the sole stdout output.
        let out: &mut dyn io::Write = if let Some(ref p) = self.output {
            file_out = std::fs::File::create(p).map_err(RsomicsError::Io)?;
            &mut file_out
        } else if self.common.json {
            sink = io::sink();
            &mut sink
        } else {
            stdout_lock = io::stdout().lock();
            &mut stdout_lock
        };
        match self.input {
            Some(ref p) => merge(p.as_path(), out),
            None => merge_stdin(out),
        }
    }
}

pub const HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Merge overlapping BED intervals (bedtools merge equivalent). Input must be sorted.",
    origin: Some(Origin {
        upstream: "bedtools",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/bioinformatics/btq033"),
    }),
    usage_lines: &["[OPTIONS] [INPUT]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: None,
                long: "input",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: false,
                default: Some("stdin"),
                description: "Input sorted BED (default: stdin)",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "output",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: false,
                default: Some("stdout"),
                description: "Output BED path (default: stdout)",
                why_default: None,
            },
            FlagSpec {
                short: Some('h'),
                long: "help",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Show this help",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Merge overlapping intervals (file must be sorted)",
            command: "rsomics-bed-sort a.bed | rsomics-bed-merge",
        },
        Example {
            description: "Merge from file",
            command: "rsomics-bed-merge sorted.bed",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use clap::CommandFactory;
    #[test]
    fn cli_definition_is_valid() {
        super::Cli::command().debug_assert();
    }
}
