use crate::turbulance_parser::{Token, TokenType, TurbulanceParseError};
use std::collections::HashMap;

/// Turbulance Language Tokenizer
/// Converts source code text into tokens for parsing
pub struct TurbulanceTokenizer {
    /// Current position in source
    position: usize,
    /// Current line number
    line: usize,
    /// Current column number
    column: usize,
    /// Source code being tokenized
    source: Vec<char>,
    /// Keyword mapping
    keywords: HashMap<String, TokenType>,
}

impl TurbulanceTokenizer {
    pub fn new(source: &str) -> Self {
        let mut keywords = HashMap::new();
        
        // Core language keywords
        keywords.insert("item".to_string(), TokenType::Item);
        keywords.insert("funxn".to_string(), TokenType::Funxn);
        keywords.insert("given".to_string(), TokenType::Given);
        keywords.insert("within".to_string(), TokenType::Within);
        keywords.insert("considering".to_string(), TokenType::Considering);
        keywords.insert("try".to_string(), TokenType::Try);
        keywords.insert("catch".to_string(), TokenType::Catch);
        keywords.insert("finally".to_string(), TokenType::Finally);
        
        // Scientific method keywords
        keywords.insert("proposition".to_string(), TokenType::Proposition);
        keywords.insert("motion".to_string(), TokenType::Motion);
        keywords.insert("evidence".to_string(), TokenType::Evidence);
        keywords.insert("metacognitive".to_string(), TokenType::Metacognitive);
        keywords.insert("goal".to_string(), TokenType::Goal);
        
        // Evidence evaluation keywords
        keywords.insert("support".to_string(), TokenType::Support);
        keywords.insert("contradict".to_string(), TokenType::Contradict);
        keywords.insert("matches".to_string(), TokenType::Matches);
        keywords.insert("contains".to_string(), TokenType::Contains);
        
        // Logical operators
        keywords.insert("and".to_string(), TokenType::And);
        keywords.insert("or".to_string(), TokenType::Or);
        keywords.insert("not".to_string(), TokenType::Not);
        
        // Literals
        keywords.insert("true".to_string(), TokenType::Boolean);
        keywords.insert("false".to_string(), TokenType::Boolean);
        keywords.insert("none".to_string(), TokenType::None);
        
        Self {
            position: 0,
            line: 1,
            column: 1,
            source: source.chars().collect(),
            keywords,
        }
    }
    
    /// Tokenize the entire source code
    pub fn tokenize(&mut self) -> Result<Vec<Token>, TurbulanceParseError> {
        let mut tokens = Vec::new();
        
        while !self.is_at_end() {
            let token = self.next_token()?;
            if token.token_type != TokenType::EOF {
                tokens.push(token);
            }
        }
        
        tokens.push(Token {
            token_type: TokenType::EOF,
            value: String::new(),
            line: self.line,
            column: self.column,
        });
        
        Ok(tokens)
    }
    
    /// Get the next token from the source
    fn next_token(&mut self) -> Result<Token, TurbulanceParseError> {
        // Skip whitespace
        self.skip_whitespace();
        
        if self.is_at_end() {
            return Ok(Token {
                token_type: TokenType::EOF,
                value: String::new(),
                line: self.line,
                column: self.column,
            });
        }
        
        let start_line = self.line;
        let start_column = self.column;
        let ch = self.advance();
        
        match ch {
            // Single-character tokens
            '(' => Ok(self.make_token(TokenType::LeftParen, "(", start_line, start_column)),
            ')' => Ok(self.make_token(TokenType::RightParen, ")", start_line, start_column)),
            '[' => Ok(self.make_token(TokenType::LeftBracket, "[", start_line, start_column)),
            ']' => Ok(self.make_token(TokenType::RightBracket, "]", start_line, start_column)),
            '{' => Ok(self.make_token(TokenType::LeftBrace, "{", start_line, start_column)),
            '}' => Ok(self.make_token(TokenType::RightBrace, "}", start_line, start_column)),
            ',' => Ok(self.make_token(TokenType::Comma, ",", start_line, start_column)),
            ';' => Ok(self.make_token(TokenType::Semicolon, ";", start_line, start_column)),
            ':' => Ok(self.make_token(TokenType::Colon, ":", start_line, start_column)),
            '.' => Ok(self.make_token(TokenType::Dot, ".", start_line, start_column)),
            '?' => Ok(self.make_token(TokenType::Question, "?", start_line, start_column)),
            '!' => {
                if self.match_char('=') {
                    Ok(self.make_token(TokenType::NotEqual, "!=", start_line, start_column))
                } else {
                    Ok(self.make_token(TokenType::Exclamation, "!", start_line, start_column))
                }
            },
            
            // Operators that might be compound
            '+' => {
                if self.match_char('=') {
                    Ok(self.make_token(TokenType::PlusAssign, "+=", start_line, start_column))
                } else {
                    Ok(self.make_token(TokenType::Plus, "+", start_line, start_column))
                }
            },
            '-' => {
                if self.match_char('=') {
                    Ok(self.make_token(TokenType::MinusAssign, "-=", start_line, start_column))
                } else if self.match_char('>') {
                    Ok(self.make_token(TokenType::Arrow, "->", start_line, start_column))
                } else {
                    Ok(self.make_token(TokenType::Minus, "-", start_line, start_column))
                }
            },
            '*' => {
                if self.match_char('=') {
                    Ok(self.make_token(TokenType::MultiplyAssign, "*=", start_line, start_column))
                } else if self.match_char('*') {
                    Ok(self.make_token(TokenType::Power, "**", start_line, start_column))
                } else {
                    Ok(self.make_token(TokenType::Multiply, "*", start_line, start_column))
                }
            },
            '/' => {
                if self.match_char('=') {
                    Ok(self.make_token(TokenType::DivideAssign, "/=", start_line, start_column))
                } else if self.match_char('/') {
                    // Single-line comment
                    self.skip_line_comment();
                    self.next_token()
                } else if self.match_char('*') {
                    // Multi-line comment
                    self.skip_block_comment()?;
                    self.next_token()
                } else {
                    Ok(self.make_token(TokenType::Divide, "/", start_line, start_column))
                }
            },
            '%' => Ok(self.make_token(TokenType::Modulo, "%", start_line, start_column)),
            '=' => {
                if self.match_char('=') {
                    Ok(self.make_token(TokenType::Equal, "==", start_line, start_column))
                } else if self.match_char('>') {
                    Ok(self.make_token(TokenType::DoubleArrow, "=>", start_line, start_column))
                } else {
                    Ok(self.make_token(TokenType::Assign, "=", start_line, start_column))
                }
            },
            '<' => {
                if self.match_char('=') {
                    Ok(self.make_token(TokenType::LessEqual, "<=", start_line, start_column))
                } else {
                    Ok(self.make_token(TokenType::Less, "<", start_line, start_column))
                }
            },
            '>' => {
                if self.match_char('=') {
                    Ok(self.make_token(TokenType::GreaterEqual, ">=", start_line, start_column))
                } else {
                    Ok(self.make_token(TokenType::Greater, ">", start_line, start_column))
                }
            },
            
            // String literals
            '"' => self.string_literal(start_line, start_column),
            '\'' => self.char_literal(start_line, start_column),
            
            // Numbers
            c if c.is_ascii_digit() => self.number_literal(c, start_line, start_column),
            
            // Identifiers and keywords
            c if c.is_alphabetic() || c == '_' => self.identifier(c, start_line, start_column),
            
            // Unexpected character
            c => Err(TurbulanceParseError {
                message: format!("Unexpected character: '{}'", c),
                line_number: Some(start_line),
                column_number: Some(start_column),
            }),
        }
    }
    
    /// Parse a string literal
    fn string_literal(&mut self, start_line: usize, start_column: usize) -> Result<Token, TurbulanceParseError> {
        let mut value = String::new();
        
        while !self.is_at_end() && self.peek() != '"' {
            if self.peek() == '\n' {
                self.line += 1;
                self.column = 1;
            }
            
            if self.peek() == '\\' {
                self.advance(); // consume backslash
                if self.is_at_end() {
                    return Err(TurbulanceParseError {
                        message: "Unterminated string literal".to_string(),
                        line_number: Some(start_line),
                        column_number: Some(start_column),
                    });
                }
                
                match self.advance() {
                    'n' => value.push('\n'),
                    't' => value.push('\t'),
                    'r' => value.push('\r'),
                    '\\' => value.push('\\'),
                    '"' => value.push('"'),
                    '\'' => value.push('\''),
                    c => {
                        value.push('\\');
                        value.push(c);
                    }
                }
            } else {
                value.push(self.advance());
            }
        }
        
        if self.is_at_end() {
            return Err(TurbulanceParseError {
                message: "Unterminated string literal".to_string(),
                line_number: Some(start_line),
                column_number: Some(start_column),
            });
        }
        
        // Consume closing quote
        self.advance();
        
        Ok(self.make_token(TokenType::String, &value, start_line, start_column))
    }
    
    /// Parse a character literal
    fn char_literal(&mut self, start_line: usize, start_column: usize) -> Result<Token, TurbulanceParseError> {
        let mut value = String::new();
        
        if self.is_at_end() || self.peek() == '\'' {
            return Err(TurbulanceParseError {
                message: "Empty character literal".to_string(),
                line_number: Some(start_line),
                column_number: Some(start_column),
            });
        }
        
        if self.peek() == '\\' {
            self.advance(); // consume backslash
            if self.is_at_end() {
                return Err(TurbulanceParseError {
                    message: "Unterminated character literal".to_string(),
                    line_number: Some(start_line),
                    column_number: Some(start_column),
                });
            }
            
            match self.advance() {
                'n' => value.push('\n'),
                't' => value.push('\t'),
                'r' => value.push('\r'),
                '\\' => value.push('\\'),
                '"' => value.push('"'),
                '\'' => value.push('\''),
                c => {
                    value.push('\\');
                    value.push(c);
                }
            }
        } else {
            value.push(self.advance());
        }
        
        if self.is_at_end() || self.peek() != '\'' {
            return Err(TurbulanceParseError {
                message: "Unterminated character literal".to_string(),
                line_number: Some(start_line),
                column_number: Some(start_column),
            });
        }
        
        // Consume closing quote
        self.advance();
        
        Ok(self.make_token(TokenType::String, &value, start_line, start_column))
    }
    
    /// Parse a number literal
    fn number_literal(&mut self, first_digit: char, start_line: usize, start_column: usize) -> Result<Token, TurbulanceParseError> {
        let mut value = String::new();
        value.push(first_digit);
        
        // Consume digits
        while !self.is_at_end() && self.peek().is_ascii_digit() {
            value.push(self.advance());
        }
        
        // Check for decimal point
        if !self.is_at_end() && self.peek() == '.' && self.peek_next().map_or(false, |c| c.is_ascii_digit()) {
            value.push(self.advance()); // consume '.'
            
            while !self.is_at_end() && self.peek().is_ascii_digit() {
                value.push(self.advance());
            }
            
            // Check for scientific notation
            if !self.is_at_end() && (self.peek() == 'e' || self.peek() == 'E') {
                value.push(self.advance());
                
                if !self.is_at_end() && (self.peek() == '+' || self.peek() == '-') {
                    value.push(self.advance());
                }
                
                if self.is_at_end() || !self.peek().is_ascii_digit() {
                    return Err(TurbulanceParseError {
                        message: "Invalid number format".to_string(),
                        line_number: Some(start_line),
                        column_number: Some(start_column),
                    });
                }
                
                while !self.is_at_end() && self.peek().is_ascii_digit() {
                    value.push(self.advance());
                }
            }
            
            Ok(self.make_token(TokenType::Float, &value, start_line, start_column))
        } else {
            // Check for scientific notation on integers
            if !self.is_at_end() && (self.peek() == 'e' || self.peek() == 'E') {
                value.push(self.advance());
                
                if !self.is_at_end() && (self.peek() == '+' || self.peek() == '-') {
                    value.push(self.advance());
                }
                
                if self.is_at_end() || !self.peek().is_ascii_digit() {
                    return Err(TurbulanceParseError {
                        message: "Invalid number format".to_string(),
                        line_number: Some(start_line),
                        column_number: Some(start_column),
                    });
                }
                
                while !self.is_at_end() && self.peek().is_ascii_digit() {
                    value.push(self.advance());
                }
                
                Ok(self.make_token(TokenType::Float, &value, start_line, start_column))
            } else {
                Ok(self.make_token(TokenType::Integer, &value, start_line, start_column))
            }
        }
    }
    
    /// Parse an identifier or keyword
    fn identifier(&mut self, first_char: char, start_line: usize, start_column: usize) -> Result<Token, TurbulanceParseError> {
        let mut value = String::new();
        value.push(first_char);
        
        while !self.is_at_end() && (self.peek().is_alphanumeric() || self.peek() == '_') {
            value.push(self.advance());
        }
        
        let token_type = self.keywords.get(&value)
            .cloned()
            .unwrap_or(TokenType::Identifier);
        
        Ok(self.make_token(token_type, &value, start_line, start_column))
    }
    
    /// Skip whitespace characters
    fn skip_whitespace(&mut self) {
        while !self.is_at_end() {
            match self.peek() {
                ' ' | '\r' | '\t' => {
                    self.advance();
                },
                '\n' => {
                    self.line += 1;
                    self.column = 1;
                    self.advance();
                },
                _ => break,
            }
        }
    }
    
    /// Skip a line comment (// to end of line)
    fn skip_line_comment(&mut self) {
        while !self.is_at_end() && self.peek() != '\n' {
            self.advance();
        }
    }
    
    /// Skip a block comment (/* to */)
    fn skip_block_comment(&mut self) -> Result<(), TurbulanceParseError> {
        let start_line = self.line;
        let start_column = self.column;
        
        while !self.is_at_end() {
            if self.peek() == '*' && self.peek_next() == Some('/') {
                self.advance(); // consume '*'
                self.advance(); // consume '/'
                return Ok(());
            }
            
            if self.peek() == '\n' {
                self.line += 1;
                self.column = 1;
            }
            
            self.advance();
        }
        
        Err(TurbulanceParseError {
            message: "Unterminated block comment".to_string(),
            line_number: Some(start_line),
            column_number: Some(start_column),
        })
    }
    
    /// Helper methods
    fn is_at_end(&self) -> bool {
        self.position >= self.source.len()
    }
    
    fn advance(&mut self) -> char {
        if !self.is_at_end() {
            self.column += 1;
            let ch = self.source[self.position];
            self.position += 1;
            ch
        } else {
            '\0'
        }
    }
    
    fn peek(&self) -> char {
        if self.is_at_end() {
            '\0'
        } else {
            self.source[self.position]
        }
    }
    
    fn peek_next(&self) -> Option<char> {
        if self.position + 1 >= self.source.len() {
            None
        } else {
            Some(self.source[self.position + 1])
        }
    }
    
    fn match_char(&mut self, expected: char) -> bool {
        if self.is_at_end() || self.peek() != expected {
            false
        } else {
            self.advance();
            true
        }
    }
    
    fn make_token(&self, token_type: TokenType, value: &str, line: usize, column: usize) -> Token {
        Token {
            token_type,
            value: value.to_string(),
            line,
            column,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_basic_tokenization() {
        let source = r#"
            item temperature = 23.5
            funxn calculate_average(numbers):
                given temperature > 30:
                    print("Hot weather")
        "#;
        
        let mut tokenizer = TurbulanceTokenizer::new(source);
        let tokens = tokenizer.tokenize().unwrap();
        
        assert_eq!(tokens[0].token_type, TokenType::Item);
        assert_eq!(tokens[1].token_type, TokenType::Identifier);
        assert_eq!(tokens[1].value, "temperature");
        assert_eq!(tokens[2].token_type, TokenType::Assign);
        assert_eq!(tokens[3].token_type, TokenType::Float);
        assert_eq!(tokens[3].value, "23.5");
    }
    
    #[test]
    fn test_proposition_tokenization() {
        let source = r#"
            proposition DrugEfficacy:
                motion ReducesSymptoms("Drug reduces symptom severity")
                
                within clinical_data:
                    given symptom_reduction > 0.5:
                        support ReducesSymptoms
        "#;
        
        let mut tokenizer = TurbulanceTokenizer::new(source);
        let tokens = tokenizer.tokenize().unwrap();
        
        // Find proposition token
        let prop_token = tokens.iter().find(|t| t.token_type == TokenType::Proposition).unwrap();
        assert_eq!(prop_token.value, "proposition");
        
        // Find motion token
        let motion_token = tokens.iter().find(|t| t.token_type == TokenType::Motion).unwrap();
        assert_eq!(motion_token.value, "motion");
        
        // Find support token
        let support_token = tokens.iter().find(|t| t.token_type == TokenType::Support).unwrap();
        assert_eq!(support_token.value, "support");
    }
    
    #[test]
    fn test_string_literals() {
        let source = r#""Hello, world!" 'x' "Line 1\nLine 2""#;
        
        let mut tokenizer = TurbulanceTokenizer::new(source);
        let tokens = tokenizer.tokenize().unwrap();
        
        assert_eq!(tokens[0].token_type, TokenType::String);
        assert_eq!(tokens[0].value, "Hello, world!");
        
        assert_eq!(tokens[1].token_type, TokenType::String);
        assert_eq!(tokens[1].value, "x");
        
        assert_eq!(tokens[2].token_type, TokenType::String);
        assert_eq!(tokens[2].value, "Line 1\nLine 2");
    }
    
    #[test]
    fn test_numbers() {
        let source = "42 3.14159 1.23e-4 2E+10";
        
        let mut tokenizer = TurbulanceTokenizer::new(source);
        let tokens = tokenizer.tokenize().unwrap();
        
        assert_eq!(tokens[0].token_type, TokenType::Integer);
        assert_eq!(tokens[0].value, "42");
        
        assert_eq!(tokens[1].token_type, TokenType::Float);
        assert_eq!(tokens[1].value, "3.14159");
        
        assert_eq!(tokens[2].token_type, TokenType::Float);
        assert_eq!(tokens[2].value, "1.23e-4");
        
        assert_eq!(tokens[3].token_type, TokenType::Float);
        assert_eq!(tokens[3].value, "2E+10");
    }
    
    #[test]
    fn test_operators() {
        let source = "+ - * / % ** == != < > <= >= += -= *= /= -> =>";
        
        let mut tokenizer = TurbulanceTokenizer::new(source);
        let tokens = tokenizer.tokenize().unwrap();
        
        let expected_types = vec![
            TokenType::Plus, TokenType::Minus, TokenType::Multiply, TokenType::Divide,
            TokenType::Modulo, TokenType::Power, TokenType::Equal, TokenType::NotEqual,
            TokenType::Less, TokenType::Greater, TokenType::LessEqual, TokenType::GreaterEqual,
            TokenType::PlusAssign, TokenType::MinusAssign, TokenType::MultiplyAssign,
            TokenType::DivideAssign, TokenType::Arrow, TokenType::DoubleArrow,
        ];
        
        for (i, expected_type) in expected_types.iter().enumerate() {
            assert_eq!(tokens[i].token_type, *expected_type);
        }
    }
    
    #[test]
    fn test_comments() {
        let source = r#"
            // This is a line comment
            item x = 42
            /* This is a
               block comment */
            item y = 24
        "#;
        
        let mut tokenizer = TurbulanceTokenizer::new(source);
        let tokens = tokenizer.tokenize().unwrap();
        
        // Should only have tokens for the actual code, not comments
        let identifiers: Vec<_> = tokens.iter()
            .filter(|t| t.token_type == TokenType::Identifier)
            .map(|t| &t.value)
            .collect();
        
        assert_eq!(identifiers, vec!["x", "y"]);
    }
}