module ripemd::ripemd160 {
    use std::u32;
    use std::u64;
    use std::vector;
    use std::string;

    /// Performs a left-rotation on a 32-bit unsigned integer.
    /// Equivalent to `(x << s) | (x >> (32 - s))`.
    fun rol_u32(x: u32, s: u8): u32 {
        ((x << s) | (x >> (32 - s)))
    }

    /// The RIPEMD-160 f function, varying per round.
    /// It returns a nonlinear function of the inputs, depending on the round number `j`.
    fun f(j: u8, x: u32, y: u32, z: u32): u32 {
        if (j < 16) {
            x ^ y ^ z
        } else if (j < 32) {
            (x & y) | (!x & z)
        } else if (j < 48) {
            (x | !y) ^ z
        } else if (j < 64) {
            (x & z) | (y & !z)
        } else {
            x ^ (y | !z)
        }
    }

    /// Returns the constant value for round `j`.
    fun K(j: u8): u32 {
        if (j < 16) 0x00000000
        else if (j < 32) 0x5A827999
        else if (j < 48) 0x6ED9EBA1
        else if (j < 64) 0x8F1BBCDC
        else 0xA953FD4E
    }

    const H0: u32 = 0x67452301;
    const H1: u32 = 0xEFCDAB89;
    const H2: u32 = 0x98BADCFE;
    const H3: u32 = 0x10325476;
    const H4: u32 = 0xC3D2E1F0;

    const R: vector<u8> = vector::from_array([
        0, 1, 2, 3, 4, 5, 6, 7,
        8, 9, 10, 11, 12, 13, 14, 15,
        7, 4, 13, 1, 10, 6, 15, 3,
        12, 0, 9, 5, 2, 14, 11, 8,
        3, 10, 14, 4, 9, 15, 8, 1,
        2, 7, 0, 6, 13, 11, 5, 12,
        1, 9, 11, 10, 0, 8, 12, 4,
        13, 3, 7, 15, 14, 5, 6, 2,
        4, 0, 5, 9, 7, 12, 2, 10,
        14, 1, 3, 8, 11, 6, 15, 13
    ]);

    const S: vector<u8> = vector::from_array([
        11, 14, 15, 12, 5, 8, 7, 9,
        11, 13, 14, 15, 6, 7, 9, 8,
        7, 6, 8, 13, 11, 9, 7, 15,
        7, 12, 15, 9, 11, 7, 13, 12,
        11, 13, 6, 7, 14, 9, 13, 15,
        14, 8, 13, 6, 5, 12, 7, 5,
        11, 12, 14, 15, 14, 15, 9, 8,
        9, 14, 5, 6, 8, 6, 5, 12,
        9, 15, 5, 11, 6, 8, 13, 12,
        5, 12, 13, 14, 11, 8, 5, 6
    ]);

    /// Computes the RIPEMD-160 hash of a given input message.
    /// This function pads the input according to the RIPEMD-160 specification
    /// and processes the data in 512-bit (64-byte) chunks.
    ///
    /// Returns a 20-byte hash (vector<u8> of length 20).
    public fun ripemd160(message: vector<u8>): vector<u8> {
        let padded = pad_input(message);

        let mut i = 0;
        let mut A = H0;
        let mut B = H1;
        let mut C = H2;
        let mut D = H3;
        let mut E = H4;

        while (i < vector::length(&padded)) {
            let block = vector::sub_range(&padded, i, i + 64);

            let mut X = vector::empty<u32>();
            let mut j = 0;
            while (j < 64) {
                let b0 = *vector::borrow(&block, j) as u32;
                let b1 = *vector::borrow(&block, j+1) as u32;
                let b2 = *vector::borrow(&block, j+2) as u32;
                let b3 = *vector::borrow(&block, j+3) as u32;
                let word = b0 | (b1 << 8) | (b2 << 16) | (b3 << 24);
                vector::push_back(&mut X, word);
                j = j + 4;
            }

            let mut AA = A;
            let mut BB = B;
            let mut CC = C;
            let mut DD = D;
            let mut EE = E;

            let mut j = 0;
            while (j < 80) {
                let r = *vector::borrow(&R, j);
                let s = *vector::borrow(&S, j);
                let T = rol_u32(AA.wrapping_add(f(j, BB, CC, DD)).wrapping_add(*vector::borrow(&X, r)).wrapping_add(K(j)), s);
                let tmp = EE;
                EE = DD;
                DD = CC;
                CC = rol_u32(BB, 10);
                BB = AA;
                AA = T.wrapping_add(tmp);
                j = j + 1;
            }

            A = A.wrapping_add(AA);
            B = B.wrapping_add(BB);
            C = C.wrapping_add(CC);
            D = D.wrapping_add(DD);
            E = E.wrapping_add(EE);

            i = i + 64;
        }

        let mut out = vector::empty<u8>();
        push_u32(&mut out, A);
        push_u32(&mut out, B);
        push_u32(&mut out, C);
        push_u32(&mut out, D);
        push_u32(&mut out, E);
        out
    }

    /// Pads the input message to a multiple of 64 bytes as required by RIPEMD-160.
    /// Appends a single 0x80 byte, then enough 0x00 bytes to make the total length
    /// congruent to 56 mod 64. Finally appends the original length (in bits) as 8 bytes, LE.
    fun pad_input(message: vector<u8>): vector<u8> {
        let len = vector::length(&message);
        let mut padded = vector::copy(&message);
        let bit_len = (len as u64) * 8;

        vector::push_back(&mut padded, 0x80);

        let mut pad_len = ((56 - ((len + 1) % 64)) + 64) % 64;
        let mut i = 0;
        while (i < pad_len) {
            vector::push_back(&mut padded, 0x00);
            i = i + 1;
        }

        push_u64_le(&mut padded, bit_len);
        padded
    }

    /// Appends a 64-bit unsigned integer to the vector in little-endian format.
    fun push_u64_le(vec: &mut vector<u8>, val: u64) {
        let mut i = 0;
        while (i < 8) {
            let byte = ((val >> (i * 8)) & 0xff) as u8;
            vector::push_back(vec, byte);
            i = i + 1;
        }
    }

    /// Appends a 32-bit unsigned integer to the vector in little-endian format.
    fun push_u32(vec: &mut vector<u8>, x: u32) {
        vector::push_back(vec, (x & 0xff) as u8);
        vector::push_back(vec, ((x >> 8) & 0xff) as u8);
        vector::push_back(vec, ((x >> 16) & 0xff) as u8);
        vector::push_back(vec, ((x >> 24) & 0xff) as u8);
    }

    /// Unit test verifying the RIPEMD-160 implementation against official test vectors.
    #[test]
    public fun test_ripemd160_vectors() {
        fun assert_eq_hex(result: vector<u8>, expected: vector<u8>) {
            let n = vector::length(&result);
            assert!(n == vector::length(&expected), 100);
            let mut i = 0;
            while (i < n) {
                assert!(*vector::borrow(&result, i) == *vector::borrow(&expected, i), 101 + i);
                i = i + 1;
            }
        }

        let msg1 = vector::empty<u8>();
        let out1 = ripemd160(msg1);
        let exp1 = vector::from_array([
            0x9c, 0x11, 0x85, 0xa5, 0xc5, 0xe9, 0xfc, 0x54,
            0x61, 0x28, 0x08, 0x97, 0x7e, 0xe8, 0xf5, 0x48,
            0xb2, 0x25, 0x8d, 0x31
        ]);
        assert_eq_hex(out1, exp1);

        let msg2 = vector::from_array([0x61, 0x62, 0x63]); // "abc"
        let out2 = ripemd160(msg2);
        let exp2 = vector::from_array([
            0x8e, 0xb2, 0x08, 0xf7, 0xe0, 0x5d, 0x98, 0x7a,
            0x9b, 0x04, 0x4a, 0x8e, 0x98, 0xc6, 0xb0, 0x87,
            0xf1, 0x5a, 0x0b, 0xfc
        ]);
        assert_eq_hex(out2, exp2);

        let msg3 = vector::from_array([
            0x6d, 0x65, 0x73, 0x73, 0x61, 0x67, 0x65, 0x20,
            0x64, 0x69, 0x67, 0x65, 0x73, 0x74
        ]); // "message digest"
        let out3 = ripemd160(msg3);
        let exp3 = vector::from_array([
            0x5d, 0x06, 0x89, 0xef, 0x49, 0xd2, 0xfa, 0xe5,
            0x72, 0xb8, 0x81, 0xb1, 0x23, 0xa8, 0x5f, 0xfa,
            0x21, 0x59, 0x5f, 0x36
        ]);
        assert_eq_hex(out3, exp3);
    }
}
