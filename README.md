Below is an example of a well‐formatted `README.md` file you can copy and paste into your repository:

---

```markdown
# Quantitative Futures Strategy

This repository contains a C++ backtesting engine for quantitative futures strategies. It downloads intraday historical data, calculates various technical indicators, and runs a backtest on commodity futures. The results (including performance metrics and trade logs) are output as CSV files for further analysis.

## Features

- **Data Fetching:**  
  Uses [libcurl](https://curl.se/libcurl/) to download historical intraday data from the Financial Modeling Prep API, and parses JSON responses with [nlohmann/json](https://github.com/nlohmann/json).

- **Indicator Calculations:**  
  Computes technical indicators such as:
  - Simple Moving Average (SMA)
  - Exponential Moving Average (EMA)
  - Standard Deviation & Z-Score
  - ATR (Average True Range) using EMA
  - Momentum signals and reaction multipliers

- **Backtesting Engine:**  
  Simulates entry and exit signals based on indicator signals and risk management rules. Tracks portfolio performance and individual trades.

- **Performance Metrics:**  
  Calculates key performance metrics including annualized return, volatility, Sharpe ratio, Sortino ratio, max drawdown, Calmar ratio, and detailed trade statistics.

- **CSV Output:**  
  Exports backtest results and trade summaries to CSV files for further review and analysis.

## Dependencies

- A C++ compiler with C++11 (or later) support
- [CMake](https://cmake.org) for build configuration
- [libcurl](https://curl.se/libcurl/) for HTTP requests
- [nlohmann/json](https://github.com/nlohmann/json) for JSON parsing

Make sure these libraries are installed on your system.

## Build Instructions

1. **Clone the repository:**

   ```bash
   git clone https://github.com/Gannon-Stoner/QuantitativeFuturesStrategy.git
   ```

2. **Navigate to the project directory:**

   ```bash
   cd QuantitativeFuturesStrategy
   ```

3. **Create a build directory and configure the project:**

   ```bash
   mkdir build && cd build
   cmake ..
   ```

4. **Build the project:**

   ```bash
   make
   ```

5. **Run the executable:**

   ```bash
   ./QuantitativeFuturesStrategy
   ```

## Usage

- **API Key:**  
  Update the `apikey` variable in `main.cpp` with your valid Financial Modeling Prep API key.

- **Commodity Symbols:**  
  The project is pre-configured with several commodity symbols (e.g., wheat, soybeans, gold, crude oil). You can modify or extend these in `main.cpp`.

- **Results:**  
  The backtest outputs CSV files containing:
  - Portfolio values and timestamps
  - Detailed trade logs (entry/exit times, prices, profit, etc.)

Performance metrics are also printed to the console upon completion.

## Code Overview

- **Data Structures:**  
  - `Bar`: Represents a single historical data point.
  - `Trade`: Stores details for individual trades.
  - `PerformanceMetrics`: Aggregates performance statistics.

- **DataFetcher Module:**  
  Downloads and parses intraday historical data using libcurl and nlohmann/json.

- **Indicator Functions:**  
  Functions for computing SMA, EMA, standard deviation, ATR, and other indicators.

- **Backtest Engine:**  
  Implements the trading logic, risk management, and portfolio simulation.

- **CSV Output:**  
  Writes backtest results and trade summaries to CSV files for analysis.

## License

This project is licensed under the MIT License.
