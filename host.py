import qsharp

from qsharp import Result
from Quantum.Bell import TestBellState

initials = (Result.Zero, Result.One)

for i in initials:
  res = TestBellState.simulate(count=1000, initial=i)
  (num_zeros, num_ones) = res
  print(f'Init:{i: <4} 0s={num_zeros: <4} 1s={num_ones: <4}')