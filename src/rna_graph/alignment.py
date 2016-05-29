class AlignmentRecord:
    NUM_FIELDS = 11

    def __init__(self, s1, e1, s2, e2, len1, len2, idy, lenr, lenq, ref,
                 query):
        self.s1, self.e1, self.s2, self.e2, self.len1, self.len2, self.idy,\
            self.lenr, self.lenq, self.ref, self.query = s1, e1, s2, e2, len1,\
            len2, idy, lenr, lenq, ref, query

    @classmethod
    def from_line(cls, line):
        fields = list(line.split())
        assert len(fields) == AlignmentRecord.NUM_FIELDS
        s1, e1, s2, e2, len1, len2 = map(int, fields[:6])
        idy = float(fields[6])
        lenr, lenq = map(int, fields[7:9])
        ref, query = fields[9:]
        return AlignmentRecord(s1, e1, s2, e2, len1, len2, idy, lenr, lenq,
                               ref, query)

    def __str__(self):
        return ' '.join(map(str, self.s1, self.e1, self.s2, self.e2, self.len1,
                            self.len2, self.idy, self.lenr, self.lenq,
                            self.ref, self.query))
