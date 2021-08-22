use ark_serialize::CanonicalSerialize;

// Add some functionality to our transcript. Everything in arkworks is CanonicalSerialize, so it'd
// be nice to just be able to pass those values directly
pub(crate) trait TranscriptProtocol {
    fn append_serializable<S: CanonicalSerialize>(&mut self, label: &'static [u8], val: &S);
}

impl TranscriptProtocol for merlin::Transcript {
    fn append_serializable<S: CanonicalSerialize>(&mut self, label: &'static [u8], val: &S) {
        let mut buf = Vec::new();
        val.serialize(&mut buf);
        self.append_message(label, &buf);
    }
}
