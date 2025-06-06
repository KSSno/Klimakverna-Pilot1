import arrow
import json
from pathlib import Path
from pydantic import BaseModel
from pydantic import RootModel
from pydantic import model_validator
from pydantic import ValidationError
from pydantic.color import Color
from typing import Literal


class Variable(BaseModel):
    name: Literal["pr", "tas"]
    description: str
    units: str

    @model_validator(mode="after")
    def check_units(self):
        if self.name == "pr":
            assert self.units == "kg m-2 s-1", "If 'variable' is pr, 'units' need to be kg m-2 s-1"
        elif self.name == "tas":
            assert self.units == "K", "If 'variable' is tas, 'units' need to be K"

        return self


class Indicator(BaseModel):
    id: int
    description: str
    units: str | None

class CalculationType(RootModel):
    root: dict[Literal["DailyTimeSeries", "30YearStatistics"], Literal["mean"]]

    @model_validator(mode="after")
    def check_keys_and_values(self):
        allowed_keys = ["DailyTimeSeries", "30YearStatistics"]
        allowed_value = "mean"

        for key, value in self.root.items():
            assert key in allowed_keys, f"Invalid calculation type '{key}'"
            assert value == allowed_value, f"Invalid statistic '{value}' for calculation type '{key}'"
        return self

class Period(BaseModel):
    description: str
    valid_period: dict[str, int]
    start: int | None = None
    end: int | None = None

    @model_validator(mode="after")
    def validate_start_end(self):
        if self.start and self.end:
            assert arrow.get(self.end).is_between(
                arrow.get(self.valid_period["start"]), arrow.get(self.valid_period["end"]), "[]"
            ) and arrow.get(self.start).is_between(
                arrow.get(self.valid_period["start"]), arrow.get(self.valid_period["end"]), "[]"
            ), "Start and end for calculation have to be in the specified period"

        if self.start and not self.end:
            assert False, "If start year is specified, end year has to be specified too"

        if self.end and not self.start:
            assert False, "If end year is specified, start year has to be specified too"

        return self

class Periods(RootModel):
    root: dict[str, Period]

    @model_validator(mode="after")
    def check_periods(self):
        for key in self.root:
            assert key in ["hist", "nf", "ff"], "Period has to be 'hist', 'nf' or 'ff'"
        return self
    
    # @model_validator(mode="after")
    # def check_valid_period(self):
    #     if self.root == "hist":
    #         assert self.root["valid_period"]["start"] == 1990 and self.root["valid_period"]["end"] == 2020, "Valid period for 'hist' is 1990-2020"
    #     if self.root == "nf":
    #         assert self.root["valid_period"]["start"] == 2041 and self.root["valid_period"]["end"] == 2070, "Valid period for 'hist' is 1990-2020"
    #     if self.root == "ff":
    #         assert self.root["valid_period"]["start"] == 2071 and self.root["valid_period"]["end"] == 2100, "Valid period for 'hist' is 1990-2020"

class Scenario(BaseModel):
    description: str
    color: str

    # @model_validator(mode="after")
    # def color_must_be_hex(self):
    #     assert self.color == Color(self.color).as_hex(), f"Color {self.color} is not a valid hex color"
    #     return self


class Scenarios(RootModel):
    root: dict[str, Scenario]


class Configuration(BaseModel):
    indicator: Indicator
    variable: Variable
    calculation_type: CalculationType
    periods: Periods
    scenarios: Scenarios
    models: list[str]
    bias_method: str
    region: dict[str, str | int]

def load_config_file(config_file_path: Path) -> dict[str, any]:
    try:
        with open(config_file_path) as f:
            json_dict = json.load(f)
        return json_dict
    except (FileNotFoundError, json.JSONDecodeError):
        raise FileNotFoundError(f"Configuration file {config_file_path} was not found or is not valid JSON")


def validate_config(config_file: str) -> dict[str, any]:
    """
    Validates the configuration file and returns a dictionary with the validated configuration

    Args:
        configuration_file: Path to the configuration file

    Returns:
        dict[str, any]: Validated configuration dictionary
    """
    
    try:
        config_content = load_config_file(Path(config_file))
        Configuration(**config_content)
        return config_content
    except (ValidationError, KeyError, TypeError) as e:
        print(f"Configuration validation failed: {e}")


# if __name__ == "__main__":
#     validated_config = validate_config("/lustre/storeC-ext/users/klimakverna/development/Klimakverna-Pilot1/config/testcase_8/config_CMIP5.json")
#     print("YEAH")