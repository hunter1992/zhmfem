use std::io::{BufWriter, Writer};

pub struct Zwriter {
    pub file_path: &str,
    pub file_writer: BufWriter,
}

impl Zwriter {
    pub fn new(file_path: &str, file_writer: BufWriter) -> Zwriter {
        Zwriter {
            file_path,
            file_writer,
        }
    }

   pub fn write_buffer(&mut self, ) 
}
