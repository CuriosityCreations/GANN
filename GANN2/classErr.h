// classErr.h defines types of errors that may be thrown and caught.
// GANN CORE
// Copyright (C) 2000-2004 Robert G. Beiko

// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// Modified from file 'classErr.h' created by Dr. Robert L. Charlebois, NeuroGadgets, Inc.

#ifndef ERROR_INCLUDED
#define ERROR_INCLUDED

#include "DCGA.h"

class Err {
protected:
	std::string message_;
public:
	explicit Err(std::string& theMessage) : message_(theMessage) { };
	Err(const Err& rhs) { message_ = rhs.message_; };
	~Err() { };
	Err& operator=(const Err& rhs) { message_ = rhs.message_; return *this; };
	const std::string& what() const { return message_; };
};

class FileErr : public Err {
public:
	explicit FileErr(std::string& theFile) : Err(theFile) { };
	~FileErr() { };
};

class FileOpenErr : public FileErr {
public:
	explicit FileOpenErr(std::string& theFile) : FileErr(theFile) { };
	~FileOpenErr() { };
};

class FileWriteErr : public FileErr {
public:
	explicit FileWriteErr(std::string& theFile) : FileErr(theFile) { };
	~FileWriteErr() { };
};

class FileReadErr : public FileErr {
public:
	explicit FileReadErr(std::string& theFile) : FileErr(theFile) { };
	~FileReadErr() { };
};

class PvmErr {
protected:
	int errorCode_;
public:
	explicit PvmErr(int theCode) : errorCode_(theCode) { };
	PvmErr(const PvmErr& rhs) { errorCode_ = rhs.errorCode_; };
	~PvmErr() { };
	PvmErr& operator=(const PvmErr& rhs) { errorCode_ = rhs.errorCode_; return *this; };
	int itsCode() const { return errorCode_; };
};

#endif
