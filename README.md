# Quantitative Futures Strategy

C++ backtester for a momentum strategy on commodity futures. Pulls intraday bars from the
Financial Modeling Prep API, computes SMA/EMA/z-score/ATR indicators, simulates long/short
entries and exits with ATR-based position sizing, and reports annualized return, volatility,
Sharpe, Sortino, Calmar, max drawdown, and per-trade statistics to the console and CSV.

## Building

```bash
mkdir build && cd build
cmake ..
make
```

Requires CMake, libcurl, and [nlohmann/json](https://github.com/nlohmann/json).

## Running

```bash
export FMP_API_KEY=your_key
./FuturesStrategyImproved
```

Symbols, timeframe, and date range are set near the bottom of `main.cpp`. Each symbol's
backtest writes a portfolio-value CSV and a trade log next to the binary, and a summary of
the performance metrics is printed to the console.

MIT license.
