"""
Unit conversion helper functions.

Provides commonly used conversion factors between different unit systems
(SI, CGS, medical units, etc.).
"""

# Length conversions
cm_to_m() = 1e-2

m_to_mm() = 1e3

# Area conversions
cm2_to_m2() = 1e-4

m2_to_cm2() = 1e4

# Volume conversions
ml_to_m3() = 1e-6

m3_to_ml() = 1e6

# Pressure conversions
Pa_to_mmHg() = 1 / 133.322

mmHg_to_Pa() = 133.322
