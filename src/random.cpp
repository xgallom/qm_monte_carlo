//
// Created by xgallom on 4/16/19.
//

#include "random.h"
#include "config.h"

Random::Random() :
		m_device(),
		m_engine(m_device())
{}

float Random::get()
{ return m_dist(m_engine); }

float Random::norm()
{ return m_norm(m_engine); }

float Random::dist()
{ return Config::Dist * get(); }
