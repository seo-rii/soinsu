# 문제 7 최종 풀이

## 정답

```text
FLAG{d1rty_b1t_l34k_c0pp3rsm1th_m33ts_str4t3gy}
```

복원한 두 소인수는 다음과 같다.

```text
p = 0xffa360d46885c534d186538179633fafc2c0548a2e24a2c1c878569f522939e38ff75ead1bcd442a834974a52e3ac66c1bc2ee00d63e2de4c436d2ca740a624699e1a1af94045c63261323cb723a7ba1b3c00bbbb6a4e534c11469e73ddbb2e50b0bc2461fcbac0726360c2c0ac9450a9a892cbf1d98ceee48827591ccc593c9
q = 0xe557fa8670389cb60c84416a65742a74fd11ed33b1631f787e92b90887b5391dacba00a386911bf8a8fbd57430b9b26e455329405ffe289e20616fe3b5562ea9b533f8f8db94bb8dcd280a6af056108e176008d3655428ad0ac6396318ba0f6efe496eac3f8585675bfed67081e0c518be5685e4daf7060abe1c58b73cc5f1e9
```

## 공격 구성

누출되지 않은 양끝 4비트 블록을 후보 번호 하나로 묶었다.

```text
low  = p[150..153]
high = p[920..923]
cid  = (high << 4) | low
```

각 후보를 고정하면 나머지 미지 비트를 두 변수로 묶어 다음 선형식을 얻는다.

```text
p = C + 2^265*x + 2^600*y
0 <= x < 2^155
0 <= y < 2^230
p | N
```

`experiments/grouped_hm_flatter.cpp`는 이 식에 대해 centered triangular Herrmann--May 격자를 만든다. 최종 스캔 파라미터는 `m=17`, `t=5`, 171차원, `x`-monic, 목표 RHF `1.15`였다. 계수는 적절한 `N` 거듭제곱으로 중심화하고, 열은 `X^i Y^j`로 스케일한 뒤 FLATTER로 감축했다.

감축 벡터에서 얻은 짧은 보조 다항식들은 FLINT로 여러 소수체 위에서 결과식을 계산했다. 오답 후보는 첫 소수 `1,000,000,007`에서 결과식들의 GCD가 1이 되어 즉시 배제됐다. 정답 후보 `cid=155`에서는 다음 값이 복원됐다.

```text
low  = 0xb
high = 0x9
```

8개 소수체에서 같은 `(x, y)` 공통근을 얻었고, CRT 모듈러스가 240비트가 되었을 때 위의 `p`를 정확히 재구성했다.

## 재현

패키지의 `experiments` 디렉터리에서 빌드한다.

```bash
g++ -O3 -std=c++17 -I/usr/local/include grouped_hm_flatter.cpp \
  -L/usr/local/lib -Wl,-rpath,/usr/local/lib \
  -lflatter -lflint -lgmpxx -lgmp -lopenblas -pthread \
  -o grouped_hm_flatter
```

정답 edge 후보는 다음 명령으로 복원할 수 있다.

```bash
OPENBLAS_NUM_THREADS=1 ./grouped_hm_flatter \
  --challenge --cid 155 --m 17 --t 5 --rhf 1.15 --lead x --centered
```

## 복호화와 검증

```python
q = N // p
phi = (p - 1) * (q - 1)
d = pow(65537, -1, phi)
m = pow(ct, d, N)
plaintext = m.to_bytes((m.bit_length() + 7) // 8, "big")
```

최종적으로 다음 조건을 모두 확인했다.

```text
p * q == N
(p & MASK) == LEAK
pow(m, 65537, N) == ct
m.hex() == 464c41477b64317274795f6231745f6c33346b5f633070703372736d3174685f6d333374735f73747234743367797d
```
