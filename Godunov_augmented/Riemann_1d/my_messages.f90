module my_messages	


real(8), parameter :: error_data=-1.0d12


contains


subroutine my_error_message

print*, 'There must be something wrong!'
pause

end subroutine my_error_message


end module my_messages
