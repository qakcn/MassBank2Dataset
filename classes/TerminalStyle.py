# This script contains classes that used to print pretty outputs.
#
# Author: qakcn
# Email: qakcn@hotmail.com
# Date: 2023-11-29

# Copyright 2023 qakcn
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


if __name__ == "__main__":
    raise SystemExit("This script is not meant to be run directly")

class TSMeta(type):
    """Terminal Style Meta Class"""
    def __getattr__(cls, name):
        def method(text):
            return cls.stylize(name, text)
        return method

class TS(metaclass=TSMeta):
    """Terminal Styles"""
    
    style = dict(
        reset = "\033[0m",
        bold = "\033[1m",
        underline = "\033[4m",
        blink = "\033[5m",
        reverse = "\033[7m",
        concealed = "\033[8m",
    )
    color = dict(
        black = "\033[30m",
        red = "\033[31m",
        green = "\033[32m",
        yellow = "\033[33m",
        blue = "\033[34m",
        magenta = "\033[35m",
        cyan = "\033[36m",
        white = "\033[37m",
    )
    bg_color = dict(
        black_bg = "\033[40m",
        red_bg = "\033[41m",
        green_bg = "\033[42m",
        yellow_bg = "\033[43m",
        blue_bg = "\033[44m",
        magenta_bg = "\033[45m",
        cyan_bg = "\033[46m",
        white_bg = "\033[47m",
    )

    printer = print

    @classmethod
    def stylize(cls, style: str, text) -> str:
        """Stylize text with style"""
        if style == cls.style["reset"]:
            return text
        if text is str and text.endswith(cls.style["reset"]):
            text = text[:-4]
        if style in cls.style:
            start = cls.style[style]
        elif style in cls.color:
            start = cls.color[style]
        elif style in cls.bg_color:
            start = cls.bg_color[style]
        text = start + str(text) + cls.style["reset"]
        return text
    
    @classmethod
    def indent(cls, text, level: int = 1, spaces: str = " "*4) -> str:
        return spaces*level + str(text)
    
    @classmethod
    def register_printer(cls, printer) -> None:
        cls.printer = printer
    
    @classmethod
    def print(cls, text, flush: bool = True, **kwargs) -> None:
        cls.printer(text, flush=flush, **kwargs)
    
    @classmethod
    def inline_print(cls, text, end: str = "", **kwargs) -> None:
        cls.print(text, end=end, **kwargs)

    @classmethod
    def indent_print(cls, text, inline = False, level: int = 1, spaces: str = " "*4, **kwargs) -> None:
        indent_text = cls.indent(cls, text, level=level, spaces=spaces)
        if inline:
            cls.inline_print(indent_text, **kwargs)
        else:
            cls.print(indent_text, **kwargs)
    
    @classmethod
    def p(cls, *args, **kwargs) -> None:
        cls.print(*args, **kwargs)

    @classmethod
    def ip(cls, *args, **kwargs) -> None:
        cls.inline_print(*args, **kwargs)
    
    @classmethod
    def idp(cls, *args, **kwargs) -> None:
        cls.indent_print(*args, **kwargs)
