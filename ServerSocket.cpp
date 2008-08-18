// Implementation of the ServerSocket class

#include "ServerSocket.h"

ServerSocket::ServerSocket ( int port )
{
  if ( ! Socket::create() )
    {
      throw "Could not create server socket." ;
    }

  if ( ! Socket::bind ( port ) )
    {
      throw "Could not bind to port." ;
    }

  if ( ! Socket::listen() )
    {
      throw "Could not listen to socket." ;
    }

}

ServerSocket::~ServerSocket()
{
}


const ServerSocket& ServerSocket::operator << ( const std::string& s ) const
{
  if ( ! Socket::send ( s ) )
    {
      throw "Could not write to socket." ;
    }

  return *this;

}


const ServerSocket& ServerSocket::operator >> ( std::string& s ) const
{

  if ( ! Socket::recv ( s ) )
    {
      throw "Could not read from socket." ;
    }

  return *this;
}

void ServerSocket::accept ( ServerSocket& sock )
{
  if ( ! Socket::accept ( sock ) )
    {
      throw "Could not accept socket." ;
    }
}

void ServerSocket::remove ( ServerSocket& sock )
{
  if ( ! Socket::remove ( sock ) )
    {
      throw "Could not remove socket." ;
    }
}
