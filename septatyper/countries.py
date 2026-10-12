"""Country name or code to ISO 3166-1 alpha-3."""
from typing import Optional

import pycountry

COUNTRY_ALIASES = {"UK": "GBR", "Russia": "RUS", "South Korea": "KOR", "Vietnam": "VNM"}


def country_to_iso3(name: Optional[str]) -> Optional[str]:
    if not name:
        return None
    if name in COUNTRY_ALIASES:
        return COUNTRY_ALIASES[name]
    try:
        return pycountry.countries.lookup(name).alpha_3
    except LookupError:
        return None
