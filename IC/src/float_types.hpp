template<typename FloatType>
struct floatinfo {
  static constexpr char const *name = "unknown";
  static constexpr int doubleprecision = 0;
};

template<>
struct floatinfo<double> {
  static constexpr char const *name = "doub";
  static constexpr int doubleprecision = 1;
};

template<>
struct floatinfo<float> {
  static constexpr char const *name = "float";
  static constexpr int doubleprecision = 0;
};
