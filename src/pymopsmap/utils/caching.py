"""Cache key computation — blake2b hash of pydantic models."""

import hashlib
import pickle
from typing import Any

from pydantic import BaseModel


def cache_key(model: BaseModel | list[BaseModel], extra: Any = None) -> str:
    """Deterministic blake2b hash of a pydantic model (or list of models)."""
    data = pickle.dumps(model, protocol=5)
    if extra is not None:
        data += pickle.dumps(extra, protocol=5)
    return hashlib.blake2b(data).hexdigest()
