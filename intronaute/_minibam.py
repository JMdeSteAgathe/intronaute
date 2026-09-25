"""Tiny embedded BAM (base64) used by the self-test.

Layout (chr1, + strand):
  exon1 [800,1200)  intron [1200,1700)  exon2 [1700,2100)
  12 correctly spliced reads (exact junction), 5 intronic reads,
  4 reads crossing the 5' boundary, 3 crossing the 3' boundary,
  + 1 duplicate / 1 MAPQ=3 / 1 secondary alignment (all must be filtered out).
A second contig named "2" (no "chr" prefix) exercises chromosome-name mapping.
"""

MINI_BAM_B64 = (
    "H4sIBAAAAAAA/wYAQkMCAGoAc3L0ZXRiYGBw8HDhDPOzMtQz4wz2t0rOzy9KycxLLEnlcggO5Az2"
    "s0rOKDLk9PGzMjIAAZioEUjIECLEBDSFFYhBKhkceJkZQAJGDAvaGBkAPZCLpmUAAAAfiwgEAAAA"
    "AAD/BgBCQwIAawPd2z1oE3EcxvGLl5c6KFgFFRRfBhV0aF7bokOhIDjo4KKDEIyNUCxaK9XgoIuD"
    "Li46uLi4COLi4uIiTetQBRcRRBFBBBFEEEEUsV7NPVBK8nvORZ7fHaQtZOmH3OX7z+9/aQed43AY"
    "BLl9B/oz0d9j0WMhPhafaxXrA8FIPgjW/POx6z8chw6OBu2YccRmFJ0wjtmMkhPGcZtRdsI4aTMq"
    "ThinbEbVCWPSZtScMKZtxqATxkWbMeSEcdlmDIszZmPGlYiRN/KnnnE4rhGHesfhuE4c6iGH4wZx"
    "qJccjlvEoZ5yOG4Th3rL4bhDHOoxh+MucajXHI77xKGeczgeEId6z+djx6VsEGQjR9jF0WwOBFuj"
    "J0a3BH9/e9cUU6UppUpTTpWmkipNNVWaWqo0g6nSDKVKM+xE83SJJte7nl4WAwk5XlYDGCBdzfYe"
    "IDXHJ9QHLwkZ6nOXhAz1sUtChvrUZSZmrMt13oe7McZPq18aUOw3FepXBhQtU6F+YXRRZJcrxqYn"
    "xRXtJYpc2P21mDhz4aw4o8uLkVnOONc84URxz7wwvLzVvjEV6uNtXBhTeat7U+rFSMhQT0ZChnoz"
    "wHhpMFolN3c4vbYZXk6qdzbDy0n1wWaoZwOMTzbDSze+2Az1XVEwvtkM9U1RMH7YDPU9UTB+2wz1"
    "LVEwwoLJUN8Rxc5uX8HY2S35ucNpFXGodxyOfuJQDzkc64lDveRwbCIO9ZTDsY041FsOxw7iUI85"
    "HLuJQ73mcAwQh3rO4agSh3rP56L/c/F//xk9Cj0cjbqDbx4BssAg6kUHJJshEPWkA7KSQdSbDshq"
    "BlGPOiBrGUS96oBsYBD1rAOymUHUuw7IdgZRDzsgOxlEvexPYsie6EefGUT1tENSpBL1tkNSoxL1"
    "uEOyl0rU6w7JCJWo5/1ZLJnJdD6TdLulr1HH9xXCjer39CX1FFPmKaXMU06Zp5IyT9WJZzb2HF3R"
    "e+bSqOvfyQjHeeJQX83AcZM41NcyWPI/DM0lv/69KYA8YhD1MwuQxwzi5dSaYxD1VTIg8wyivkgG"
    "5DmDeJmBvWAQLzOwVwziZQb2lkG8zMDeM4iXGdjH0Pxs7+B+FUg+U4l62yH5SiXqcYfkO5Wo1x2S"
    "X1Sinvc/gF3IHEhTAAAfiwgEAAAAAAD/BgBCQwIAGwADAAAAAAAAAAAA"
)

MINI_BAI_B64 = (
    "QkFJAQIAAAACAAAASRIAAAEAAAAAAGsAAAAAAHg0awAAAAAASpIAAAIAAAAAAGsAAAAAAHg0awAA"
    "AAAAQwAAAAAAAAAAAAAAAAAAAAEAAAAAAGsAAAAAAAIAAABJEgAAAQAAAHg0awAAAAAAAADXAwAA"
    "AABKkgAAAgAAAHg0awAAAAAAAADXAwAAAAAnAAAAAAAAAAAAAAAAAAAAAQAAAHg0awAAAAAAAAAA"
    "AAAAAAA="
)
