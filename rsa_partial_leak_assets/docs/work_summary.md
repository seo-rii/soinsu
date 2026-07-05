# 2026 암호분석경진대회 7번 RSA 부분 비트 유출 문제 작업 정리

작성 목적: 외부 실행 환경에서 작업을 이어가기 위한 재현 가능한 코드와 현재까지의 분석 내용을 정리한다.

## 1. 문제 요약

문제는 2048비트 RSA 모듈러스 `N = p*q`와 공개지수 `e = 65537`, 암호문 `ct`, 그리고 1024비트 소수 `p`의 일부 비트 누출 정보 `MASK`, `p & MASK`가 주어진 상황이다. `MASK`의 1비트 위치가 `p`의 알려진 비트이고, 목표는 `N`을 인수분해하여 textbook RSA `ct = m^e mod N`을 복호화한 뒤 평문 바이트열(big-endian)을 출력하는 것이다.

## 2. 비트 구조 분석

`MASK`를 LSB 기준으로 분석하면 다음과 같다.

```text
p의 알려진 비트 수: 672 / 1024
p의 미지 비트 수:   352 / 1024

미지 블록:
bits 150..153    길이 4
bits 265..348    길이 84
bits 362..419    길이 58
bits 600..668    길이 69
bits 682..768    길이 87
bits 784..829    길이 46
bits 920..923    길이 4
```

즉 알려진 비트 비율은 65.625%이다. 미지 비트는 무작위로 완전히 흩어진 것이 아니라 7개의 연속 블록에 모여 있다. 양끝의 4비트 블록 두 개를 브루트포스하면 총 256개 후보가 생기며, 나머지는 두 개의 큰 변수로 묶을 수 있다.

```text
low  = p[150..153]   4비트 후보
high = p[920..923]   4비트 후보
x    = p[265..419]   155비트 그룹
y    = p[600..829]   230비트 그룹
p    = p0 + 2^265*x + 2^600*y
```

이때 `x`와 `y` 안에는 원래 알려진 gap 비트도 포함되어 있으므로, 후보 루트가 나온 뒤 최종 검증은 반드시 다음 조건으로 수행한다.

```python
N % p == 0
(p & MASK) == LEAK
```

## 3. 이론적 방향

이 문제는 Coppersmith 계열의 부분 키 노출 공격이다. 단일 연속 미지 블록이라면 고전적인 univariate Coppersmith로 접근할 수 있지만, 여기서는 알려진 비트가 여러 위치에 있으므로 Herrmann-May식 “unknown divisor modulo linear equation” 모델이 더 자연스럽다.

선형식은 다음과 같다.

```text
f(x, y) = p0 + 2^265*x + 2^600*y
f(x, y) == 0 mod p,  p | N
```

`p` 자체는 모르지만 `N`이 `p`의 배수이므로, cuso의 `modulus='p'`, `modulus_multiple=N` 인터페이스나 Herrmann-May 격자를 사용할 수 있다. 두 변수 구성에서 루트 bound의 곱은 `2^(155+230) = 2^385`, 즉 `N^(385/2048) ≈ N^0.188` 수준이다. Herrmann-May의 bivariate unknown-divisor bound와 비교하면 이론적으로 여유가 있어 보이므로, 현재 최우선 경로는 2변수 cuso 실행이다.

## 4. 지금까지 시도한 구현

### 4.1 마스크 파싱 및 구조 확인

`constants.py`, `solve_bits.py`, `inspect_problem.py` 등으로 PDF의 `MASK`, `p & MASK`, `N`, `ct`를 파싱하고 위의 7개 미지 블록 구조를 확인했다.

### 4.2 SAT 접근

초기에는 `p`의 알려진 비트를 Boolean 상수로 고정하고, `p*q=N`의 전체 곱셈 회로를 CNF/XOR 형태로 만들었다. 전체 인코딩은 너무 커서 일반 솔버로 밀기 어려웠다. 이후 양끝 4비트 후보를 고정하여 `p mod 2^265`, `q mod 2^265`, `q`의 상위 prefix를 유도하는 축소 SAT도 만들었다. 대략 후보 하나당 100만 변수급 문제가 되어 병렬 솔버 없이는 시간이 컸다.

관련 파일:

```text
sat_factor.py
sat_candidate_opt.py
sat_mul.py
sat_mul_cnf.py
run_sat_one.py
run_sat_cnf_one.py
z3_factor.py
```

### 4.3 5변수 Herrmann-May

양끝 4비트 블록을 브루트포스한 뒤 나머지 5개 미지 블록을 각각 변수로 두는 격자도 구성했다.

```text
p = p0 + 2^265*x1 + 2^362*x2 + 2^600*x3 + 2^682*x4 + 2^784*x5
```

`m=3, t=2`는 56차원 정도로 빠르게 실행되었지만 충분히 강한 보조 다항식이 나오지 않았다. `m=4, t=2` 이상은 현재 컨테이너의 fpylll 환경에서 오래 걸렸다.

관련 파일:

```text
hm_multi.py
hm_multivar.py
hm5_lattice_check.py
basis_m3.txt
basis_m4.txt
```

### 4.4 2변수 Herrmann-May 및 루트 복원

가장 유망한 방향이다. `x = p[265..419]`, `y = p[600..829]`로 묶고, 256개 후보 각각에 대해 bivariate lattice를 생성했다. LLL에서 얻은 보조 다항식의 공통근은 작은 소수체 위에서 roots를 찾고 CRT로 합치는 방식과 Groebner/resultant 계열 방식 모두 실험했다.

현재 컨테이너에는 Sage/cuso/flatter가 없어 cuso 경로는 실행하지 못했다. pure fpylll fallback으로 `m=7,t=2`, `m=8,t=3` 등을 일부 시험했지만, 아직 `p`를 찾지 못했다.

관련 파일:

```text
solve7_main.py
solve7_compact_previous.py
hm_bivar.py
hm_core.py
hm_scan_gb.py
root_recover.py
root_groebner.py
```

## 5. 패키지의 권장 실행법

패키지 루트에서 다음을 실행한다.

```bash
sage -python src/solve7_main.py --mode analyze
sage -python src/solve7_main.py --mode cuso --a 0 --b 256 | tee logs/cuso_full.log
```

후보를 나누어 실행하려면 다음처럼 한다.

```bash
bash scripts/run_cuso_range.sh 0 64
bash scripts/run_cuso_range.sh 64 128
bash scripts/run_cuso_range.sh 128 192
bash scripts/run_cuso_range.sh 192 256
```

8개 병렬 작업은 다음처럼 실행할 수 있다.

```bash
bash scripts/run_cuso_8way.sh
```

pure Python fallback은 다음과 같다.

```bash
python3 -m pip install -r env/requirements-local.txt
python3 src/solve7_main.py --mode local --m 8 --t 3 --lead y --a 0 --b 1
```

단, fallback은 검증 및 실험용이다. 실제 해결은 Sage+cuso+flatter 환경이 더 적합하다.

## 6. 결과가 나온 뒤 복호화

`solve7_main.py`는 후보 `p`가 검증되면 자동으로 복호화한다. 외부 도구로 `p`만 얻은 경우에는 다음을 실행한다.

```bash
python3 src/decrypt_with_factor.py <p-hex-or-decimal>
```

내부적으로는 다음 표준 RSA 복호화를 수행한다.

```python
q = N // p
phi = (p - 1) * (q - 1)
d = pow(65537, -1, phi)
m = pow(ct, d, N)
plaintext = m.to_bytes((m.bit_length() + 7) // 8, 'big')
```

## 7. 패키지 구성

```text
problem/
  problem_7_partial_factorization.pdf
src/
  solve7_main.py
  inspect_problem.py
  decrypt_with_factor.py
scripts/
  run_cuso_range.sh
  run_cuso_8way.sh
  run_local_range.sh
env/
  setup_external.md
  requirements-local.txt
  Dockerfile.local
docs/
  work_summary.docx
  work_summary.md
experiments/
  중간 실험 코드와 로그
logs/
  실행 로그 저장 위치
```

## 8. 남은 작업 우선순위

1. SageMath 9.8+, cuso, flatter, msolve 환경에서 `--mode cuso` 전체 후보를 실행한다.
2. 실패하면 cuso 파라미터 또는 변수 그룹을 조정한다. 특히 `x-leading`, `y-leading`, 변수 범위에 알려진 gap 비트를 명시하는 변형을 시험한다.
3. 그래도 실패하면 SAT+Coppersmith hybrid로 넘어간다. SAT는 임의 위치의 알려진 비트를 잘 활용하고, Coppersmith는 대수 구조를 활용하므로 두 방식을 결합하는 전략이 적합하다.
4. 백업으로 kionactf/coppersmith의 multivariate linear Herrmann-May 구현을 사용해 같은 식을 재구성한다.

## 9. 참고 자료

- Mathias Herrmann, Alexander May, “Solving Linear Equations Modulo Divisors: On Factoring Given Any Bits”, ASIACRYPT 2008. https://link.springer.com/chapter/10.1007/978-3-540-89255-7_25
- Yameen Ajani, Curtis Bright, “SAT and Lattice Reduction for Integer Factorization”, ISSAC 2024. https://arxiv.org/abs/2406.20071
- cuso EUROCRYPT 2025 artifact, “Solving Multivariate Coppersmith Problems”. https://artifacts.iacr.org/eurocrypt/2025/a13/readme.html
- kionactf/coppersmith. https://github.com/kionactf/coppersmith
- jvdsn/crypto-attacks. https://github.com/jvdsn/crypto-attacks

## 10. 현재 상태에 대한 정직한 메모

이 패키지는 “완성된 풀이”가 아니라 “외부 강한 환경에서 계속 실행하기 위한 작업 에셋”이다. 현재 컨테이너에서는 Sage/cuso/flatter 부재와 fpylll 성능 한계 때문에 최종 `p`, `q`, 평문을 얻지 못했다. 그래도 문제 상수, 비트 구조, 공격식, 후보 스캔 방식, 복호화 루틴은 모두 재현 가능하게 정리되어 있다.
